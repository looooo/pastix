/**
 *
 * @file core_zgelrops.c
 *
 * PaStiX low-rank kernel routines
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Gregoire Pichon
 * @date 2016-03-23
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <cblas.h>
#include <lapacke.h>
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "z_nan_check.h"
#include "kernels_trace.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t zzero =  0.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Allocate a low-rank matrix.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          Number of rows of the matrix A.
 *
 * @param[in] N
 *          Number of columns of the matrix A.
 *
 * @param[in] rkmax
 *         @arg -1: the matrix is allocated tight to its rank.
 *         @arg >0: the matrix is allocated to the minimum of rkmax and its maximum rank.
 *
 * @param[out] A
 *          The allocated low-rank matrix
 *
 *******************************************************************************/
void
core_zlralloc( pastix_int_t      M,
               pastix_int_t      N,
               pastix_int_t      rkmax,
               pastix_lrblock_t *A )
{
    pastix_complex64_t *u, *v;

    if ( rkmax == -1 ) {
        u = malloc( M * N * sizeof(pastix_complex64_t) );
        memset(u, 0, M * N * sizeof(pastix_complex64_t) );
        A->rk = -1;
        A->rkmax = M;
        A->u = u;
        A->v = NULL;
    }
    else {
        pastix_int_t rk = pastix_imin( M, N );
        rkmax = pastix_imin( rkmax, rk );

#if defined(PASTIX_DEBUG_LR)
        u = malloc( M * rkmax * sizeof(pastix_complex64_t) );
        v = malloc( N * rkmax * sizeof(pastix_complex64_t) );

        /* To avoid uninitialised values in valgrind. Lapacke doc (xgesvd) is not correct */
        memset(u, 0, M * rkmax * sizeof(pastix_complex64_t));
        memset(v, 0, N * rkmax * sizeof(pastix_complex64_t));
#else
        u = malloc( (M+N) * rkmax * sizeof(pastix_complex64_t));

        /* To avoid uninitialised values in valgrind. Lapacke doc (xgesvd) is not correct */
        memset(u, 0, (M+N) * rkmax * sizeof(pastix_complex64_t));

        v = u + M * rkmax;
#endif

        A->rk = 0;
        A->rkmax = rkmax;
        A->u = u;
        A->v = v;
    }
}

/**
 *******************************************************************************
 *
 * @brief Free a low-rank matrix.
 *
 *******************************************************************************
 *
 * @param[inout] A
 *          The low-rank matrix that will be desallocated.
 *
 *******************************************************************************/
void
core_zlrfree( pastix_lrblock_t *A )
{
    if ( A->rk == -1 ) {
        free(A->u);
        A->u = NULL;
    }
    else {
        free(A->u);
#if defined(PASTIX_DEBUG_LR)
        free(A->v);
#endif
        A->u = NULL;
        A->v = NULL;
    }
    A->rk = 0;
    A->rkmax = 0;
}

/**
 *******************************************************************************
 *
 * @brief Resize a low-rank matrix
 *
 *******************************************************************************
 *
 * @param[in] copy
 *          Enable/disable the copy of the data from A->u and A->v into the new
 *          low-rank representation.
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix A.
 *
 * @param[inout] A
 *          The low-rank representation of the matrix. At exit, this structure
 *          is modified with the new low-rank representation of A, is the rank
 *          is small enough
 *
 * @param[in] newrk
 *          The new rank of the matrix A.
 *
 * @param[in] newrkmax
 *          The new maximum rank of the matrix A. Useful if the low-rank
 *          structure was allocated with more data than the rank.
 *
 *******************************************************************************
 *
 * @return  The new rank of A
 *
 *******************************************************************************/
int
core_zlrsze( int copy, pastix_int_t M, pastix_int_t N,
             pastix_lrblock_t *A,
             int newrk, int newrkmax )
{
    pastix_int_t minmn = pastix_imin( M, N );

    newrkmax = (newrkmax == -1) ? newrk : newrkmax;

    /**
     * It is not interesting to compress, so we alloc space to store the full matrix
     */
    if ( (newrk * 2) > minmn )
    {
        A->u = realloc( A->u, M * N * sizeof(pastix_complex64_t) );
#if defined(PASTIX_DEBUG_LR)
        free(A->v);
#endif
        A->v = NULL;
        A->rk = -1;
        A->rkmax = M;
        return -1;
    }
    /**
     * The rank is null, we free everything
     */
    else if (newrkmax == 0)
    {
        /**
         * The rank is nul, we free everything
         */
        free(A->u);
#if defined(PASTIX_DEBUG_LR)
        free(A->v);
#endif
        A->u = NULL;
        A->v = NULL;
        A->rkmax = newrkmax;
        A->rk = newrk;
    }
    /**
     * The rank is non null, we allocate the correct amount of space, and
     * compress the stored information if necessary
     */
    else {
        pastix_complex64_t *u, *v;
        int ret;

        if ( newrkmax != A->rkmax ) {
#if defined(PASTIX_DEBUG_LR)
            u = malloc( M * newrkmax * sizeof(pastix_complex64_t) );
            v = malloc( N * newrkmax * sizeof(pastix_complex64_t) );
#else
            u = malloc( (M+N) * newrkmax * sizeof(pastix_complex64_t) );
            v = u + M * newrkmax;
#endif
            if ( copy ) {
                ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, newrk,
                                           A->u, M, u, M );
                assert(ret == 0);
                ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', newrk, N,
                                           A->v, A->rkmax, v, newrkmax );
                assert(ret == 0);
            }
            free(A->u);
#if defined(PASTIX_DEBUG_LR)
            free(A->v);
#endif
            A->u = u;
            A->v = v;
        }
        A->rk = newrk;
        A->rkmax = newrkmax;

        (void)ret;
    }
    assert( A->rk <= A->rkmax);
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Convert a low rank matrix into a dense matrix.
 *
 * Convert a low-rank matrix of size m-by-n into a full rank matrix.
 * A = op( u * v^t ) with op(A) = A or A^t
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          @arg PastixNoTrans: returns A = u * v^t
 *          @arg PastixTrans: returns A = v * u^t
 *
 * @param[in] m
 *          Number of rows of the low-rank matrix Alr.
 *
 * @param[in] n
 *          Number of columns of the low-rank matrix Alr.
 *
 * @param[in] Alr
 *          The low rank matrix to be converted into a full rank matrix
 *
 * @param[inout] A
 *          The matrix of dimension lda-by-k in which to store the uncompressed
 *          version of Alr. k = n if trans == PastixNoTrans, m otherwise.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1, m) if trans ==
 *          PastixNoTrans, lda >= max(1,n) otherwise.
 *
 *******************************************************************************
 *
 * @retval  0  in case of success.
 * @retval  -i if the ith parameter is incorrect.
 *
 *******************************************************************************/
int
core_zlr2ge( pastix_trans_t trans, pastix_int_t m, pastix_int_t n,
             const pastix_lrblock_t *Alr,
             pastix_complex64_t *A, pastix_int_t lda )
{
    int ret = 0;

#if !defined(NDEBUG)
    if ( m < 0 ) {
        return -1;
    }
    if ( n < 0 ) {
        return -2;
    }
    if (Alr == NULL || Alr->rk > Alr->rkmax) {
        return -3;
    }
    if ( (trans == PastixNoTrans && lda < m) ||
         (trans != PastixNoTrans && lda < n) )
    {
        return -5;
    }
    if ( Alr->rk == -1 ) {
        if (Alr->u == NULL || Alr->v != NULL || (Alr->rkmax < m))
        {
            return -6;
        }
    }
    else if ( Alr->rk != 0){
        if (Alr->u == NULL || Alr->v == NULL) {
            return -6;
        }
    }
#endif

    if ( trans == PastixNoTrans ) {
        if ( Alr->rk == -1 ) {
            ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                                       Alr->u, Alr->rkmax, A, lda );
            assert( ret == 0 );
        }
        else if ( Alr->rk == 0 ) {
            ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', m, n,
                                       0.0, 0.0, A, lda );
            assert( ret == 0 );
        }
        else {
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        m, n, Alr->rk,
                        CBLAS_SADDR(zone),  Alr->u, m,
                                            Alr->v, Alr->rkmax,
                        CBLAS_SADDR(zzero), A, lda);
        }
    }
    else {
        if ( Alr->rk == -1 ) {
            core_zgetro( m, n, Alr->u, Alr->rkmax, A, lda );
        }
        else if ( Alr->rk == 0 ) {
            ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', n, m,
                                       0.0, 0.0, A, lda );
            assert( ret == 0 );
        }
        else {
            cblas_zgemm(CblasColMajor, CblasTrans, CblasTrans,
                        n, m, Alr->rk,
                        CBLAS_SADDR(zone),  Alr->v, Alr->rkmax,
                                            Alr->u, m,
                        CBLAS_SADDR(zzero), A,      lda);
        }
    }

    return ret;
}

/**
 *******************************************************************************
 *
 * @brief Add the dense matrix A to the low-rank matrix B.
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 * @param[in] alpha
 *          alpha * A is add to B
 *
 * @param[in] M1
 *          The number of rows of the matrix A.
 *
 * @param[in] N1
 *          The number of columns of the matrix A.
 *
 * @param[in] A
 *          The matrix A (dense)
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *
 * @param[in] M2
 *          The number of rows of the matrix B.
 *
 * @param[in] N2
 *          The number of columns of the matrix B.
 *
 * @param[inout] B
 *          The low-rank representation of the matrix B.
 *
 * @param[in] offx
 *          The horizontal offset of A with respect to B.
 *
 * @param[in] offy
 *          The vertical offset of A with respect to B.
 *
 *******************************************************************************
 *
 * @return The new rank of B or -1 if ranks are too large for recompression
 *
 *******************************************************************************/
int
core_zgradd( const pastix_lr_t *lowrank, pastix_complex64_t alpha,
             pastix_int_t M1, pastix_int_t N1, const pastix_complex64_t *A, pastix_int_t lda,
             pastix_int_t M2, pastix_int_t N2, pastix_lrblock_t *B,
             pastix_int_t offx, pastix_int_t offy)
{
    pastix_int_t rmax = pastix_imin( M2, N2 );
    pastix_int_t rank, ldub;
    double tol = lowrank->tolerance;

    assert( B->rk <= B->rkmax);

    if ( B->rk == -1 ) {
        pastix_complex64_t *tmp = B->u;
        ldub = B->rkmax;
        tmp += ldub * offy + offx;
        core_zgeadd( PastixNoTrans, M1, N1,
                     alpha, A,   lda,
                       1.0, tmp, ldub );
        return 0;
    }

    ldub = M2;

    /**
     * The rank is too big, we need to uncompress/compress B
     */
    if ( ((B->rk + M1) > rmax) &&
         ((B->rk + N1) > rmax) )
    {
        pastix_complex64_t *work = malloc( M2 * N2 * sizeof(pastix_complex64_t) );
        assert(B->rk > 0);

        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     M2, N2, B->rk,
                     CBLAS_SADDR(zone),  B->u, ldub,
                                         B->v, B->rkmax,
                     CBLAS_SADDR(zzero), work, M2 );

        core_zgeadd( PastixNoTrans, M1, N1,
                     alpha, A, lda,
                     1.0,    work + M2 * offy + offx, M2 );

        core_zlrfree(B);
        lowrank->core_ge2lr( tol, M2, N2, work, M2, B );
        rank = B->rk;
        free(work);
    }
    /**
     * We consider the A matrix as Id * A or A *Id
     */
    else {
        pastix_lrblock_t lrA;
        lrA.rk = -1;
        lrA.rkmax = lda;
        lrA.u = (void*)A;
        lrA.v = NULL;
        rank = lowrank->core_rradd( lowrank, PastixNoTrans, &alpha,
                                    M1, N1, &lrA,
                                    M2, N2, B,
                                    offx, offy );
    }

    assert( B->rk <= B->rkmax);
    return rank;
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_lr_null
 *
 * @brief Compute the product of two possible low-rank matrices and returns the
 * result in AB.
 *
 *******************************************************************************
 * @param[in] transA
 *         @arg PastixNoTrans: No transpose, op( A ) = A;
 *         @arg PastixTrans:   Transpose, op( A ) = A';
 *
 * @param[in] transB
 *         @arg PastixNoTrans: No transpose, op( B ) = B;
 *         @arg PastixTrans:   Transpose, op( B ) = B';
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix B.
 *
 * @param[in] K
 *          The number of columns of the matrix A and the number of rows of the
 *          matrix B.
 *
 * @param[in] A
 *          The low-rank representation of the matrix A.
 *
 * @param[in] B
 *          The low-rank representation of the matrix B.
 *
 * @param[out] AB
 *          The low-rank representation of the matrix AB.
 *
 * @param[inout] work
 *          The workspace used to store temporary data
 *
 * @param[in] ldwork
 *          Dimension of the workspace work. It must be at least:
 *          - if A and B are low-rank:
 *             ldwork >= min( A->rk * ( N + B->rk ), B->rk * ( M + A->rk ) )
 *          - if A is low-rank and B full-rank
 *             ldwork >= A->rk * N
 *          - if B is low-rank and A full-rank
 *             ldwork >= B->rk * M
 *          - if A and B are full-rank
 *             ldwork >= M * N
 *
 *******************************************************************************
 *
 * @return The way the product AB is stored: AB or op(AB).
 *
 *******************************************************************************/
pastix_trans_t
core_zlrm2( pastix_trans_t transA, pastix_trans_t transB,
            pastix_int_t M, pastix_int_t N, pastix_int_t K,
            const pastix_lrblock_t *A,
            const pastix_lrblock_t *B,
            pastix_lrblock_t *AB,
            pastix_complex64_t *work,
            pastix_int_t ldwork )
{
    pastix_int_t ldau, ldav, ldbu, ldbv;
    pastix_trans_t transV = PastixNoTrans;

    assert( A->rk  <= A->rkmax);
    assert( B->rk  <= B->rkmax);
    assert(transA == PastixNoTrans);
    assert(transB != PastixNoTrans);

    /* Quick return if multiplication by 0 */
    if ( A->rk == 0 || B->rk == 0 ) {
        AB->rk = 0;
        AB->rkmax = 0;
        AB->u = NULL;
        AB->v = NULL;
        return transV;
    }

    ldau = (A->rk == -1) ? A->rkmax : M;
    ldav = A->rkmax;
    ldbu = (B->rk == -1) ? B->rkmax : N;
    ldbv = B->rkmax;

    if ( A->rk != -1 ) {
        /**
         * A and B are both low rank
         */
        if ( B->rk != -1 ) {
            errorPrint("core_zlrm2: we should not end up here\n");
            assert(0);
        }
        /**
         * A is low rank and not B
         */
        else {
            /**
             * Let's compute A * B' = Au Av^h B' by computing only Av^h * B'
             *   ABu = Au
             *   ABv = Av^h B'
             */
            assert( (A->rk * N) <= ldwork );
            AB->rk = A->rk;
            AB->rkmax = A->rk;
            AB->u = A->u;
            AB->v = work;

            kernel_trace_start_lvl2( PastixKernelLvl2_LR_GEMM_PRODUCT );
            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                         A->rk, N, K,
                         CBLAS_SADDR(zone),  A->v,  ldav,
                                             B->u,  ldbu,
                         CBLAS_SADDR(zzero), AB->v, AB->rkmax );
            kernel_trace_stop_lvl2( FLOPS_ZGEMM( A->rk, N, K ) );
        }
    }
    else {
        /**
         * B is low rank and not A
         */
        if ( B->rk != -1 ) {
            /**
             * Let's compute A * B' = A * (Bu Bv^h)' by computing only A * Bv^h'
             *   ABu = A * Bv^h'
             *   ABv = Bu'
             */
            assert( (B->rk * M) <= ldwork );
            AB->rk = B->rk;
            AB->rkmax = B->rk;
            AB->u = work;
            AB->v = B->u;

            kernel_trace_start_lvl2( PastixKernelLvl2_LR_GEMM_PRODUCT );
            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                         M, B->rk, K,
                         CBLAS_SADDR(zone),  A->u,  ldau,
                                             B->v,  ldbv,
                         CBLAS_SADDR(zzero), AB->u, M );
            kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, B->rk, K ) );

            transV = transB;
        }
        /**
         * A and B are both full rank
         */
        else {
            /**
             * A and B are both full
             *  Let's compute the product to add the full matrix
             * TODO: return low rank matrix:
             *             AB.u = A, AB.v = B' when K is small
             */
            /* if ( 2*K < pastix_imin( M, N ) ) { */
            /*     AB->rk = K; */
            /*     AB->rkmax = B->rkmax; */
            /*     AB->u = A->u; */
            /*     AB->v = B->u; */
            /*     transV = transB; */
            /* } */
            /* else { */
                assert( (M * N) <= ldwork );
                AB->rk = -1;
                AB->rkmax = M;
                AB->u = work;
                AB->v = NULL;

                kernel_trace_start_lvl2( PastixKernelLvl2_FR_GEMM );
                cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                             M, N, K,
                             CBLAS_SADDR(zone),  A->u, ldau,
                                                 B->u, ldbu,
                             CBLAS_SADDR(zzero), work, M );
                kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, N, K ) );
            /* } */
        }
    }
    assert( AB->rk <= AB->rkmax);

    /* Not used for now */
    (void)transA; (void)transB;
    (void)ldwork;

    return transV;
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_lr_null
 *
 * @brief Compute the product of two low-rank matrices (rank != -1) and returns
 * the result in AB
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 * @param[in] transA
 *         @arg PastixNoTrans: No transpose, op( A ) = A;
 *         @arg PastixTrans:   Transpose, op( A ) = A';
 *
 * @param[in] transB
 *         @arg PastixNoTrans: No transpose, op( B ) = B;
 *         @arg PastixTrans:   Transpose, op( B ) = B';
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix B.
 *
 * @param[in] K
 *          The number of columns of the matrix A and the number of rows of the
 *          matrix B.
 *
 * @param[in] A
 *          The low-rank representation of the matrix A.
 *
 * @param[in] B
 *          The low-rank representation of the matrix B.
 *
 * @param[out] AB
 *          The low-rank representation of the matrix op(A)op(B).
 *
 *******************************************************************************
 *
 * @return The way the product AB is stored: AB or op(AB).
 *
 *******************************************************************************/
pastix_trans_t
core_zlrm3( const pastix_lr_t *lowrank,
            pastix_trans_t transA, pastix_trans_t transB,
            pastix_int_t M, pastix_int_t N, pastix_int_t K,
            const pastix_lrblock_t *A,
            const pastix_lrblock_t *B,
            pastix_lrblock_t *AB )
{
    pastix_int_t ldau, ldav, ldbu, ldbv;
    int transV = PastixNoTrans;
    pastix_complex64_t *work2;
    pastix_lrblock_t rArB;
    pastix_int_t flops = 0;

    assert( A->rk  <= A->rkmax);
    assert( B->rk  <= B->rkmax);
    assert(transA == PastixNoTrans);
    assert(transB != PastixNoTrans);

    /* Quick return if multiplication by 0 */
    if ( A->rk == 0 || B->rk == 0 ) {
        AB->rk = 0;
        AB->rkmax = 0;
        AB->u = NULL;
        AB->v = NULL;
        return transV;
    }

    ldau = (A->rk == -1) ? A->rkmax : M;
    ldav = A->rkmax;
    ldbu = (B->rk == -1) ? B->rkmax : N;
    ldbv = B->rkmax;

    work2 = malloc( A->rk * B->rk * sizeof(pastix_complex64_t));


    /**
     * Let's compute A * B' = Au Av^h (Bu Bv^h)' with the smallest ws
     */
    kernel_trace_start_lvl2( PastixKernelLvl2_LR_GEMM_PRODUCT );
    cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                 A->rk, B->rk, K,
                 CBLAS_SADDR(zone),  A->v, ldav,
                                     B->v, ldbv,
                 CBLAS_SADDR(zzero), work2, A->rk );

    kernel_trace_stop_lvl2( FLOPS_ZGEMM( A->rk, B->rk, K ) );

    /**
     * Try to compress (Av^h Bv^h')
     */
    lowrank->core_ge2lr( lowrank->tolerance, A->rk, B->rk, work2, A->rk, &rArB );

    kernel_trace_start_lvl2( PastixKernelLvl2_LR_GEMM_PRODUCT );

    /**
     * The rank of AB is not smaller than min(rankA, rankB)
     */
    if (rArB.rk == -1){
        if ( A->rk < B->rk ) {
            /**
             *    ABu = Au
             *    ABv = (Av^h Bv^h') * Bu'
             */
            pastix_complex64_t *work = malloc( A->rk * N * sizeof(pastix_complex64_t));

            //assert( (A->rk * ( N + B->rk )) <= lwork );
            AB->rk = A->rk;
            AB->rkmax = A->rk;
            AB->u = A->u;
            AB->v = work;

            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                         A->rk, N, B->rk,
                         CBLAS_SADDR(zone),  work2,  A->rk,
                         B->u,  ldbu,
                         CBLAS_SADDR(zzero), AB->v, AB->rkmax );
            flops = FLOPS_ZGEMM( A->rk, N, B->rk );
        }
        else {
            /**
             *    ABu = Au * (Av^h Bv^h')
             *    ABv = Bu'
             */
            pastix_complex64_t *work = malloc( B->rk * M * sizeof(pastix_complex64_t));

            //assert( (B->rk * ( M + A->rk )) <= lwork );
            AB->rk = B->rk;
            AB->rkmax = B->rk;
            AB->u = work;
            AB->v = B->u;

            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         M, B->rk, A->rk,
                         CBLAS_SADDR(zone),  A->u,  ldau,
                         work2,  A->rk,
                         CBLAS_SADDR(zzero), AB->u, M );
            flops = FLOPS_ZGEMM( M, B->rk, A->rk );

            transV = transB;

            /* free(work); */
        }
    }
    else if (rArB.rk == 0){
        AB->rk    = 0;
        AB->rkmax = 0;
        AB->u = NULL;
        AB->v = NULL;
    }
    /**
     * The rank of AB is smaller than min(rankA, rankB)
     */
    else{
        pastix_complex64_t *work = malloc( (M + N) * rArB.rk * sizeof(pastix_complex64_t));

        AB->rk    = rArB.rk;
        AB->rkmax = rArB.rk;
        AB->u = work;
        AB->v = work + M * rArB.rk;

        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     M, rArB.rk, A->rk,
                     CBLAS_SADDR(zone),  A->u,   ldau,
                                         rArB.u, A->rk,
                     CBLAS_SADDR(zzero), AB->u,  M );

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     rArB.rk, N, B->rk,
                     CBLAS_SADDR(zone),  rArB.v, rArB.rkmax,
                                         B->u, ldbu,
                     CBLAS_SADDR(zzero), AB->v, rArB.rk );

        flops = FLOPS_ZGEMM( M, rArB.rk, A->rk ) + FLOPS_ZGEMM( rArB.rk, N, B->rk );

        /* free(work); */
    }
    core_zlrfree(&rArB);
    free(work2);

    /* Not used for now */
    (void)transA;

    kernel_trace_stop_lvl2( flops );

    return transV;
}

/**
 *******************************************************************************
 *
 * @brief Compute A * B + C with three low-rank matrices
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 * @param[in] transA
 *         @arg PastixNoTrans: No transpose, op( A ) = A;
 *         @arg PastixTrans:   Transpose, op( A ) = A';
 *
 * @param[in] transB
 *         @arg PastixNoTrans: No transpose, op( B ) = B;
 *         @arg PastixTrans:   Transpose, op( B ) = B';
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix B.
 *
 * @param[in] K
 *          The number of columns of the matrix A and the number of rows of the
 *          matrix B.
 *
 * @param[in] Cm
 *          The number of rows of the matrix C.
 *
 * @param[in] Cn
 *          The number of columns of the matrix C.
 *
 * @param[in] offx
 *          The horizontal offset of A with respect to B.
 *
 * @param[in] offy
 *          The vertical offset of A with respect to B.
 *
 * @param[in] alpha
 *          The multiplier parameter: C = beta * C + alpha * AB
 *
 * @param[in] A
 *          The low-rank representation of the matrix A.
 *
 * @param[in] B
 *          The low-rank representation of the matrix B.
 *
 * @param[in] beta
 *          The addition parameter: C = beta * C + alpha * AB
 *
 * @param[inout] C
 *          The low-rank representation of the matrix C
 *
 * @param[in] work
 *          The workspace used to store temporary data
 *
 * @param[in] ldwork
 *          The leading dimension of work
 *
 * @param[in] fcblk
 *          The facing supernode
 *
 *******************************************************************************/
void
core_zlrmm( const pastix_lr_t *lowrank,
            pastix_trans_t transA, pastix_trans_t transB,
            pastix_int_t M, pastix_int_t N, pastix_int_t K,
            pastix_int_t Cm, pastix_int_t Cn,
            pastix_int_t offx, pastix_int_t offy,
            pastix_complex64_t alpha, const pastix_lrblock_t *A,
                                      const pastix_lrblock_t *B,
            pastix_complex64_t beta,        pastix_lrblock_t *C,
            pastix_complex64_t *work, pastix_int_t ldwork,
            SolverCblk *fcblk )
{
    pastix_complex64_t *tmp = NULL;
    pastix_lrblock_t AB;
    pastix_int_t ldabu, ldabv, ldcu, ldcv;
    pastix_int_t required = 0;
    int transV;
    int allocated = 0;
    double tol = lowrank->tolerance;

    assert(transA == PastixNoTrans);
    assert(transB != PastixNoTrans);
    assert( A->rk <= A->rkmax);
    assert( B->rk <= B->rkmax);
    assert( C->rk <= C->rkmax);

    /* Quick return if multiplication by 0 */
    if ( A->rk == 0 || B->rk == 0 ) {
        return;
    }

    if ( A->rk != -1 ) {
        if ( B->rk != -1 ) {
            required = 0;
            /* required = pastix_imin( A->rk * ( N + B->rk ), */
            /*                         B->rk * ( M + A->rk ) ); */
        }
        else {
            required = A->rk * N;
        }
    }
    else {
        if ( B->rk != -1 ) {
            required = B->rk * M;
        }
        else {
            required = M * N;
        }
    }

    if ( required <= ldwork ) {
        tmp = work;
    }
    else {
        tmp = malloc( required * sizeof(pastix_complex64_t));
        allocated = 1;
    }

    if (A->rk != -1 && B->rk != -1){
        /* Note that in this case tmp is not used anymore */
        /* For instance, AB.rk != -1 */
        transV = core_zlrm3( lowrank, transA, transB, M, N, K,
                             A, B, &AB );

        if (AB.rk == 0){
            if ( allocated ) {
                free(tmp);
            }
            return;
        }
    }
    else{
        transV = core_zlrm2( transA, transB, M, N, K,
                             A, B, &AB, tmp, required );
    }

    pastix_cblk_lock( fcblk );
    ldabu = (AB.rk == -1) ? AB.rkmax : M;
    ldabv = (transV == PastixNoTrans) ? AB.rkmax : N;
    ldcu = (C->rk == -1) ? C->rkmax : Cm;
    ldcv = C->rkmax;

    /**
     * The destination matrix is full rank
     */
    if (C->rk == -1) {
        pastix_complex64_t *Cptr = C->u;
        Cptr += ldcu * offy + offx;

        if ( AB.rk == -1 ) {
            core_zgeadd( PastixNoTrans, M, N,
                         alpha, AB.u, ldabu,
                         beta,  Cptr, ldcu );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_FR_GEMM );
            cblas_zgemm( CblasColMajor, CblasNoTrans, transV,
                         M, N, AB.rk,
                         CBLAS_SADDR(alpha), AB.u, ldabu,
                                             AB.v, ldabv,
                         CBLAS_SADDR(beta),  Cptr, ldcu );
            kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, N, AB.rk ) );
        }
    }
     /**
     * The destination matrix is low rank
     */
    else {
        if ( AB.rk == -1 ) {
            assert(beta == 1.0);
            core_zgradd( lowrank, alpha,
                         M, N, tmp, AB.rkmax,
                         Cm, Cn, C,
                         offx, offy );
        }
        else {
            if ( AB.rk + C->rk > pastix_imin(M, N) ) {
                pastix_complex64_t *work = malloc( Cm * Cn * sizeof(pastix_complex64_t) );

                /* Do not uncompress a null LR structure */
                if (C->rk > 0){
                    /* Uncompress C */
                    kernel_trace_start_lvl2( PastixKernelLvl2_LR_UNCOMPRESS );
                    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                                 Cm, Cn, C->rk,
                                 CBLAS_SADDR(beta),  C->u, ldcu,
                                                     C->v, ldcv,
                                 CBLAS_SADDR(zzero), work, Cm );
                    kernel_trace_stop_lvl2( FLOPS_ZGEMM( Cm, Cn, C->rk ) );
                }
                else{
                    memset(work, 0, Cm * Cn * sizeof(pastix_complex64_t) );
                }

                /* Add A*B */
                kernel_trace_start_lvl2( PastixKernelLvl2_FR_GEMM );
                cblas_zgemm( CblasColMajor, CblasNoTrans, transV,
                             M, N, AB.rk,
                             CBLAS_SADDR(alpha), AB.u, ldabu,
                                                 AB.v, ldabv,
                             CBLAS_SADDR(zone), work + Cm * offy + offx, Cm );
                kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, N, AB.rk ) );

                core_zlrfree(C);
                lowrank->core_ge2lr( tol, Cm, Cn, work, Cm, C );
                free(work);
            }
            else {
                /* Need to handle correctly this case */
                lowrank->core_rradd( lowrank, transV, &alpha,
                                     M, N, &AB,
                                     Cm, Cn, C,
                                     offx, offy );
            }
        }
    }
    pastix_cblk_unlock( fcblk );

    /* Free memory from zlrm3 */
    if (A->rk != -1 && B->rk != -1){
        if (AB.rk == A->rk || AB.rk == B->rk){
            if ( A->rk < B->rk ) {
                free(AB.v);
            }
            else{
                free(AB.u);
            }
        }
        else if (AB.rk == 0 || A->rk == 0 || B->rk == 0){

        }
        else {
            free(AB.u);
        }
    }

    if ( allocated ) {
        free(tmp);
    }

    assert( C->rk <= C->rkmax);
}

/**
 *******************************************************************************
 *
 * @brief Compute A * B + C with A, and B low-rank matrices, and C full rank
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 * @param[in] transA
 *         @arg PastixNoTrans: No transpose, op( A ) = A;
 *         @arg PastixTrans:   Transpose, op( A ) = A';
 *
 * @param[in] transB
 *         @arg PastixNoTrans: No transpose, op( B ) = B;
 *         @arg PastixTrans:   Transpose, op( B ) = B';
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix B.
 *
 * @param[in] K
 *          The number of columns of the matrix A and the number of rows of the
 *          matrix B.
 *
 * @param[in] alpha
 *          The multiplier parameter: C = beta * C + alpha * AB
 *
 * @param[in] A
 *          The low-rank representation of the matrix A.
 *
 * @param[in] B
 *          The low-rank representation of the matrix B.
 *
 * @param[in] beta
 *          The addition parameter: C = beta * C + alpha * AB
 *
 * @param[inout] C
 *          The matrix C (full).
 *
 * @param[in] ldc
 *          The leading dimension of the matrix C.
 *
 * @param[in] work
 *          The workspace used to store temporary data
 *
 * @param[in] ldwork
 *          The leading dimension of work
 *
 * @param[in] fcblk
 *          The facing supernode
 *
 *******************************************************************************/
void
core_zlrmge( const pastix_lr_t *lowrank,
             pastix_trans_t transA, pastix_trans_t transB,
             pastix_int_t M, pastix_int_t N, pastix_int_t K,
             pastix_complex64_t alpha, const pastix_lrblock_t *A,
                                       const pastix_lrblock_t *B,
             pastix_complex64_t beta, pastix_complex64_t *C, int ldc,
             pastix_complex64_t *work, pastix_int_t ldwork,
             SolverCblk *fcblk )
{
    pastix_lrblock_t lrC;

    lrC.rk = -1;
    lrC.rkmax = ldc;
    lrC.u = C;
    lrC.v = NULL;

    core_zlrmm( lowrank, transA, transB, M, N, K,
                M, N, 0, 0,
                alpha, A, B, beta, &lrC,
                work, ldwork,
                fcblk );
}


/**
 *******************************************************************************
 *
 * @brief Copy a small low-rank structure into a large one
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 * @param[in] transA
 *         @arg PastixNoTrans: No transpose, op( A ) = A;
 *         @arg PastixTrans:   Transpose, op( A ) = A';
 *
 * @param[in] transB
 *         @arg PastixNoTrans: No transpose, op( B ) = B;
 *         @arg PastixTrans:   Transpose, op( B ) = B';
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix B.
 *
 * @param[in] K
 *          The number of columns of the matrix A and the number of rows of the
 *          matrix B.
 *
 * @param[in] alpha
 *          The multiplier parameter: C = beta * C + alpha * AB
 *
 * @param[in] A
 *          The low-rank representation of the matrix A.
 *
 * @param[in] B
 *          The low-rank representation of the matrix B.
 *
 * @param[in] beta
 *          The addition parameter: C = beta * C + alpha * AB
 *
 * @param[inout] C
 *          The matrix C (full).
 *
 * @param[in] ldc
 *          The leading dimension of the matrix C.
 *
 * @param[in] work
 *          The workspace used to store temporary data
 *
 * @param[in] ldwork
 *          The leading dimension of work
 *
 * @param[in] fcblk
 *          The facing supernode
 *
 *******************************************************************************/
void
core_zlrcpy( const pastix_lr_t *lowrank,
             pastix_trans_t transA, pastix_complex64_t alpha,
             pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
             pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
             pastix_int_t offx, pastix_int_t offy,
             pastix_complex64_t *work, pastix_int_t ldwork )
{
    pastix_complex64_t *u, *v;
    pastix_int_t ldau, ldav;
    int ret;

    assert( B->rk == 0 );
    assert( (M1 + offx) <= M2 );
    assert( (N1 + offy) <= N2 );

    ldau = (A->rk == -1) ? A->rkmax : M1;
    ldav = (transA == PastixNoTrans) ? A->rkmax : N1;

    if ( A->rk == -1 ) {
        /*
         * TODO: This case can be improved by compressing A, and then
         * copying it into B, however the criteria to keep A compressed or
         * not must be based on B dimension, and not on A ones
         */
        if ( ldwork < M2 *N2 ) {
            MALLOC_INTERN( work, M2 * N2, pastix_complex64_t );
            ldwork = -1;
        }

        if ( M1 != M2 || N1 != N2 ) {
            LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M2, N2,
                                 0.0, 0.0, work, M2 );
        }
        ret = core_zgeadd( PastixNoTrans, M1, N1,
                           alpha, A->u, ldau,
                           0.0, work + M2 * offy + offx, M2 );
        assert(ret == 0);

        lowrank->core_ge2lr( lowrank->tolerance, M2, N2, work, M2, B );

        if (ldwork == -1) {
            memFree_null(work);
        }
    }
    else {
        core_zlralloc( M2, N2, A->rkmax, B );
        u = B->u;
        v = B->v;
        B->rk = A->rk;

        if ( M1 != M2 ) {
            LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M2, B->rk,
                                 0.0, 0.0, u, M2 );
        }
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M1, A->rk,
                                   A->u, ldau,
                                   u + offx, M2 );
        assert(ret == 0);

        if ( N1 != N2 ) {
            LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', B->rk, N2,
                                 0.0, 0.0, v, B->rkmax );
        }
        ret = core_zgeadd( transA, A->rk, N1,
                           alpha, A->v, ldav,
                           0.0, v + B->rkmax * offy, B->rkmax );
        assert(ret == 0);
    }
    assert( B->rk <= B->rkmax);
    (void) ret;
}


/**
 *******************************************************************************
 *
 * @brief Concatenate left parts of two low-rank matrices
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          alpha * A is add to B
 *
 * @param[in] M1
 *          The number of rows of the matrix A.
 *
 * @param[in] N1
 *          The number of columns of the matrix A.
 *
 * @param[in] A
 *          The low-rank representation of the matrix A.
 *
 * @param[in] M2
 *          The number of rows of the matrix B.
 *
 * @param[in] B
 *          The low-rank representation of the matrix B.
 *
 * @param[in] offx
 *          The horizontal offset of A with respect to B.
 *
 * @param[inout] u1u2
 *          The workspace where matrices are concatenated
 *
 *******************************************************************************/
void
core_zlrconcatenate_u( pastix_complex64_t alpha,
                       pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                       pastix_int_t M2,                        pastix_lrblock_t *B,
                       pastix_int_t offx,
                       pastix_complex64_t *u1u2 )
{
    pastix_complex64_t *tmp;
    pastix_int_t i, ret, rank;
    pastix_int_t ldau, ldbu;

    rank = (A->rk == -1) ? pastix_imin(M1, N1) : A->rk;
    rank += B->rk;

    ldau = (A->rk == -1) ? A->rkmax : M1;
    ldbu = M2;

    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M2, B->rk,
                               B->u, ldbu, u1u2, M2 );
    assert(ret == 0);

    tmp = u1u2 + B->rk * M2;
    if ( A->rk == -1 ) {
        /*
         * A is full of rank M1, so A will be integrated into v1v2
         */
        if ( M1 < N1 ) {
            if ( M1 != M2 ) {
                /* Set to 0 */
                memset(tmp, 0, M2 * M1 * sizeof(pastix_complex64_t));

                /* Set diagonal */
                tmp += offx;
                for (i=0; i<M1; i++, tmp += M2+1) {
                    *tmp = 1.0;
                }
            }
            else {
                assert( offx == 0 );
                ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M2, M1,
                                           0.0, 1.0, tmp, M2 );
                assert( ret == 0 );
            }
        }
        else {
            /*
             * A is full of rank N1, so A is integrated into u1u2
             */
            if ( M1 != M2 ) {
                memset(tmp, 0, M2 * N1 * sizeof(pastix_complex64_t));
            }
            ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M1, N1,
                                       A->u, ldau, tmp + offx, M2 );
            assert(ret == 0);
        }
    }
    /*
     * A is low rank of rank A->rk
     */
    else {
        if ( M1 != M2 ) {
            memset(tmp, 0, M2 * A->rk * sizeof(pastix_complex64_t));
        }
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M1, A->rk,
                                   A->u, ldau, tmp + offx, M2 );
        assert(ret == 0);
    }
    (void) ret;
    (void) alpha;
}


/**
 *******************************************************************************
 *
 * @brief Concatenate right parts of two low-rank matrices
 *
 *******************************************************************************
 *
 * @param[in] transA1
 *         @arg PastixNoTrans:  No transpose, op( A ) = A;
 *         @arg PastixTrans:  Transpose, op( A ) = A';
 *
 * @param[in] alpha
 *          alpha * A is add to B
 *
 * @param[in] M1
 *          The number of rows of the matrix A.
 *
 * @param[in] N1
 *          The number of columns of the matrix A.
 *
 * @param[in] A
 *          The low-rank representation of the matrix A.
 *
 * @param[in] N2
 *          The number of columns of the matrix B.
 *
 * @param[in] B
 *          The low-rank representation of the matrix B.
 *
 * @param[in] offy
 *          The vertical offset of A with respect to B.
 *
 * @param[inout] v1v2
 *          The workspace where matrices are concatenated
 *
 *******************************************************************************/
void
core_zlrconcatenate_v( pastix_trans_t transA1, pastix_complex64_t alpha,
                       pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                                        pastix_int_t N2,       pastix_lrblock_t *B,
                       pastix_int_t offy,
                       pastix_complex64_t *v1v2 )
{
    pastix_complex64_t *tmp;
    pastix_int_t i, ret, rank;
    pastix_int_t ldau, ldav, ldbv;

    rank = (A->rk == -1) ? pastix_imin(M1, N1) : A->rk;
    rank += B->rk;

    ldau = (A->rk == -1) ? A->rkmax : M1;
    ldav = (transA1 == PastixNoTrans) ? A->rkmax : N1;
    ldbv = B->rkmax;

    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', B->rk, N2,
                               B->v, ldbv, v1v2, rank );
    assert(ret == 0);

    tmp = v1v2 + B->rk;
    if ( A->rk == -1 ) {
        assert( transA1 == PastixNoTrans );
        /*
         * A is full of rank M1, so it is integrated into v1v2
         */
        if ( M1 < N1 ) {
            if (N1 != N2) {
                ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M1, N2,
                                           0.0, 0.0, tmp, rank );
                assert( ret == 0 );
            }
            core_zgeadd( PastixNoTrans, M1, N1,
                         alpha, A->u, ldau,
                         0.0, tmp + offy * rank, rank );
        }
        /*
         * A is full of rank N1, so it has been integrated into u1u2
         */
        else {
            if (N1 != N2) {
                /* Set to 0 */
                ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', N1, N2,
                                           0.0, 0.0, tmp, rank );
                assert(ret == 0);

                /* Set diagonal */
                tmp += offy * rank;
                for (i=0; i<N1; i++, tmp += rank+1) {
                    *tmp = alpha;
                }
            }
            else {
                assert( offy == 0 );
                ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', N1, N2,
                                           0.0, alpha, tmp + offy * rank, rank );
                assert( ret == 0 );
            }
        }
    }
    /*
     * A is low rank of rank A->rk
     */
    else {
        if (N1 != N2) {
            ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', A->rk, N2,
                                       0.0, 0.0, tmp, rank );
            assert(ret == 0);
        }
        core_zgeadd( transA1, A->rk, N1,
                     alpha, A->v,              ldav,
                       0.0, tmp + offy * rank, rank );
    }
    (void) ret;
    (void) alpha;
}
