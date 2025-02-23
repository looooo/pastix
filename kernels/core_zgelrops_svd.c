/**
 *
 * @file core_zgelrops_svd.c
 *
 * PaStiX low-rank kernel routines using SVD based on Lapack ZGESVD.
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Gregoire Pichon
 * @author Esragul Korkmaz
 * @author Mathieu Faverge
 * @author Nolan Bredel
 * @date 2024-07-05
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <cblas.h>
#include <lapacke.h>
#include "flops.h"
#include "kernels_trace.h"
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"
#include "z_nan_check.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define PASTIX_SVD_2NORM 1
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#if !defined(PASTIX_SVD_2NORM)
#include "common/frobeniusupdate.h"

/**
 *******************************************************************************
 *
 * @ingroup kernel_lr_svd_null
 *
 * @brief Compute the frobenius norm of a vector
 *
 * This routine is inspired from LAPACK zlassq function, and allows to
 * accumulate the contribution backward for better accuracy as opposed to dnrm2
 * which allows only positive increment.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of elemnts in the vector
 *
 * @param[in] x
 *          The vector of size n * incx
 *
 * @param[in] incx
 *          The increment between two elments in the vector x.
 *
 *******************************************************************************
 *
 * @return  The frobenius norm of the vector x.
 *
 *******************************************************************************/
static inline double
core_dlassq( int           n,
             const double *x,
             int           incx )
{
    double scale = 1.;
    double sumsq = 0.;
    int i;

    for( i=0; i<n; i++, x+=incx ) {
        frobenius_update( 1, &scale, &sumsq, x );
    }

    return scale * sqrt( sumsq );
}
#endif /* !defined(PASTIX_SVD_2NORM) */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t zzero =  0.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Convert a full rank matrix in a low rank matrix, using SVD.
 *
 *******************************************************************************
 *
 * @param[in] use_reltol
 *          TODO
 *
 * @param[in] tol
 *          The tolerance used as a criterion to eliminate information from the
 *          full rank matrix.
 *          If tol < 0, then we compress up to rklimit. So if rklimit is set to
 *          min(m,n), and tol < 0., we get a full representation of the matrix
 *          under the form U * V^t.
 *
 * @param[in] rklimit
 *          The maximum rank to store the matrix in low-rank format. If
 *          -1, set to min(M, N) / PASTIX_LR_MINRATIO.
 *
 * @param[in] m
 *          Number of rows of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] n
 *          Number of columns of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] Avoid
 *          The matrix of dimension lda-by-n that needs to be compressed
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1, m)
 *
 * @param[out] Alr
 *          The low rank matrix structure that will store the low rank
 *          representation of A
 *
 *******************************************************************************
 *
 * @return  TODO
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zge2lr_svd( int               use_reltol,
                 pastix_fixdbl_t   tol,
                 pastix_int_t      rklimit,
                 pastix_int_t      m,
                 pastix_int_t      n,
                 const void       *Avoid,
                 pastix_int_t      lda,
                 pastix_lrblock_t *Alr )
{
    const pastix_complex64_t *A = (const pastix_complex64_t*)Avoid;
    pastix_fixdbl_t flops = 0.0;
    pastix_complex64_t *u, *v, *zwork, *Acpy, ws;
    double             *rwork, *s;
    pastix_int_t        i, ret, ldu, ldv;
    pastix_int_t        minMN, imax;
    pastix_int_t        lwork = -1;
    pastix_int_t        zsize, rsize;
    double              norm;

#if !defined(NDEBUG)
    if ( m < 0 ) {
        return -2;
    }
    if ( n < 0 ) {
        return -3;
    }
    if ( lda < m ) {
        return -5;
    }
#endif

    norm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                A, lda, NULL );

    /* Quick return on norm */
    if ( (norm == 0.) && (tol >= 0.) ) {
        core_zlralloc( m, n, 0, Alr );
        return 0. ;
    }

    rklimit = ( rklimit < 0 ) ? core_get_rklimit( m, n ) : rklimit;
    if ( tol < 0. ) {
        tol = -1.;
    }
    else if ( use_reltol ) {
        tol = tol * norm;
    }

    /* Quick return on max rank */
    minMN = pastix_imin(m, n);
    rklimit = pastix_imin( minMN, rklimit );

    /**
     * If maximum rank is 0, then either the matrix norm is below the tolerance,
     * and we can return a null rank matrix, or it is not and we need to return
     * a full rank matrix.
     */
    if ( rklimit == 0 ) {
        if ( (tol < 0.) || (norm < tol) ) {
            core_zlralloc( m, n, 0, Alr );
            return 0.;
        }

        /* Return full rank */
        core_zlralloc( m, n, -1, Alr );
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                                   A, lda, Alr->u, Alr->rkmax );
        assert(ret == 0);
        return 0.;
    }

    /*
     * Allocate a temporary Low rank matrix to store the full U and V
     */
    core_zlralloc( m, n, pastix_imin( m, n ), Alr );
    u = Alr->u;
    v = Alr->v;
    ldu = m;
    ldv = Alr->rkmax;

    /*
     * Query the workspace needed for the gesvd
     */
#if defined(PASTIX_DEBUG_LR_NANCHECK)
    ws = minMN;
#else
    {
        /*
         * rworkfix is a fix to an internal mkl bug
         * see: https://community.intel.com/t5/Intel-oneAPI-Math-Kernel-Library/MKL-2021-4-CoreDump-with-LAPACKE-cgesvd-work-or-LAPACKE-zgesvd/m-p/1341228
         */
        double rwork;
        ret = MYLAPACKE_zgesvd_work( LAPACK_COL_MAJOR, 'S', 'S',
                                     m, n, NULL, m,
                                     NULL, NULL, ldu, NULL, ldv,
                                     &ws, lwork, &rwork );
        (void)rwork;
    }
#endif
    lwork = ws;
    zsize = ws;
    zsize += m * n; /* Copy of the matrix A */

    rsize = minMN;
#if defined(PRECISION_z) || defined(PRECISION_c)
    rsize += 5 * minMN;
#endif

    zwork = malloc( zsize * sizeof(pastix_complex64_t) + rsize * sizeof(double) );
    rwork = (double*)(zwork + zsize);

    Acpy = zwork + lwork;
    s    = rwork;

    /*
     * Backup the original matrix before to overwrite it with the SVD
     */
    ret = LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', m, n,
                              A, lda, Acpy, m );
    assert( ret == 0 );

    ret = MYLAPACKE_zgesvd_work( LAPACK_COL_MAJOR, 'S', 'S',
                                 m, n, Acpy, m,
                                 s, u, ldu, v, ldv,
                                 zwork, lwork, rwork + minMN );
    if ( ret != 0 ) {
        pastix_print_error( "SVD Failed\n" );
    }

    /* Let's stop i before going too far */
    imax = pastix_imin( minMN, rklimit+1 );
    for (i=0; i<imax; i++, v+=1) {
        double frob_norm;

        /*
         * There are two different stopping criteria for SVD to decide the
         * compression rank:
         *	1) The 2-norm:
         *         Compare the singular values to the threshold
         *      2) The Frobenius norm:
         *         Compare the Frobenius norm of the trailing singular values to
         *         the threshold. Note that we use a reverse accumulation of the
         *         singular values to avoid accuracy issues.
         */
#if defined(PASTIX_SVD_2NORM)
        frob_norm = s[i];
#else
        frob_norm = core_dlassq( minMN-i, s + minMN - 1, -1 );
#endif

        if (frob_norm < tol) {
            break;
        }
        cblas_zdscal(n, s[i], v, ldv);
    }

    /*
     * Resize the space used by the low rank matrix
     */
    core_zlrsze( 1, m, n, Alr, i, -1, rklimit );

    /*
     * It was not interesting to compress, so we restore the dense version in Alr
     */
    if (Alr->rk == -1) {
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                                   A, lda, Alr->u, Alr->rkmax );
        assert(ret == 0);
    }

    (void)ret;
    memFree_null(zwork);
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Add two LR structures A=(-u1) v1^T and B=u2 v2^T into u2 v2^T
 *
 *    u2v2^T - u1v1^T = (u2 u1) (v2 v1)^T
 *    Compute QR decomposition of (u2 u1) = Q1 R1
 *    Compute QR decomposition of (v2 v1) = Q2 R2
 *    Compute SVD of R1 R2^T = u sigma v^T
 *    Final solution is (Q1 u sigma^[1/2]) (Q2 v sigma^[1/2])^T
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 * @param[in] transA1
 *         @arg PastixNoTrans: No transpose, op( A ) = A;
 *         @arg PastixTrans:   Transpose, op( A ) = A';
 *
 * @param[in] alphaptr
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
 * @param[in] N2
 *          The number of columns of the matrix B.
 *
 * @param[in] B
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
 * @return  The new rank of u2 v2^T or -1 if ranks are too large for recompression
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zrradd_svd( const pastix_lr_t      *lowrank,
                 pastix_trans_t          transA1,
                 const void             *alphaptr,
                 pastix_int_t            M1,
                 pastix_int_t            N1,
                 const pastix_lrblock_t *A,
                 pastix_int_t            M2,
                 pastix_int_t            N2,
                 pastix_lrblock_t       *B,
                 pastix_int_t            offx,
                 pastix_int_t            offy)
{
    pastix_int_t rank, M, N, minU, minV;
    pastix_int_t i, ret, lwork, new_rank;
    pastix_int_t ldau, ldav, ldbu, ldbv;
    pastix_complex64_t *u1u2, *v1v2, *R, *u, *v;
    pastix_complex64_t *tmp, *zbuf, *tauU, *tauV;
    pastix_complex64_t  querysize;
    pastix_complex64_t  alpha = *((pastix_complex64_t*)alphaptr);
    double *s;
    size_t wzsize, wdsize;
    double tol = lowrank->tolerance;
    pastix_fixdbl_t flops, total_flops = 0.0;

    rank = (A->rk == -1) ? pastix_imin(M1, N1) : A->rk;
    rank += B->rk;
    M = pastix_imax(M2, M1);
    N = pastix_imax(N2, N1);
    minU = pastix_imin(M, rank);
    minV = pastix_imin(N, rank);

    assert(M2 == M && N2 == N);
    assert(B->rk != -1);

    assert( A->rk <= A->rkmax);
    assert( B->rk <= B->rkmax);

    if ( ((M1 + offx) > M2) ||
         ((N1 + offy) > N2) )
    {
        pastix_print_error( "Dimensions are not correct" );
        assert(0 /* Incorrect dimensions */);
        return total_flops;
    }

    /*
     * A is rank null, nothing to do
     */
    if (A->rk == 0) {
        return total_flops;
    }

    ldau = (A->rk == -1) ? A->rkmax : M1;
    ldav = (transA1 == PastixNoTrans) ? A->rkmax : N1;
    ldbu = M;
    ldbv = B->rkmax;

    /*
     * Let's handle case where B is a null matrix
     *   B = alpha A
     */
    if (B->rk == 0) {
        core_zlrcpy( lowrank, transA1, alpha,
                     M1, N1, A, M2, N2, B,
                     offx, offy );
        return total_flops;
    }

    /*
     * The rank is too big, let's try to compress
     */
    if ( rank > pastix_imin( M, N ) ) {
        assert(0);
    }

    /*
     * Let's compute the size of the workspace
     */
    /* u1u2 and v1v2 */
    wzsize = (M+N) * rank;
    /* tauU and tauV */
    wzsize += minU + minV;

    /* Storage of R, u and v */
    wzsize += 3 * rank * rank;

    /* QR/LQ workspace */
    lwork = pastix_imax( M, N ) * 32;

    /* Workspace needed for the gesvd */
#if defined(PASTIX_DEBUG_LR_NANCHECK)
    querysize = rank;
#else
    {
        /*
         * rworkfix is a fix to an internal mkl bug
         * see: https://community.intel.com/t5/Intel-oneAPI-Math-Kernel-Library/MKL-2021-4-CoreDump-with-LAPACKE-cgesvd-work-or-LAPACKE-zgesvd/m-p/1341228
         */
        double rwork;
        ret = MYLAPACKE_zgesvd_work( LAPACK_COL_MAJOR, 'S', 'S',
                                     rank, rank, NULL, rank,
                                     NULL, NULL, rank, NULL, rank,
                                     &querysize, -1, &rwork );
        (void)rwork;
    }
#endif
    lwork = pastix_imax( lwork, querysize );
    wzsize += lwork;

    wdsize = rank;
#if defined(PRECISION_z) || defined(PRECISION_c)
    wdsize += 5 * rank;
#endif

    zbuf = malloc( wzsize * sizeof(pastix_complex64_t) + wdsize * sizeof(double) );
    s    = (double*)(zbuf + wzsize);

    u1u2 = zbuf + lwork;
    tauU = u1u2 + M * rank;
    v1v2 = tauU + minU;
    tauV = v1v2 + N * rank;
    R    = tauV + minV;

    /*
     * Concatenate U2 and U1 in u1u2
     *  [ u2  0  ]
     *  [ u2  u1 ]
     *  [ u2  0  ]
     */
    core_zlrconcatenate_u( alpha,
                           M1, N1, A,
                           M2,     B,
                           offx, u1u2 );

    /*
     * Perform QR factorization on u1u2 = (Q1 R1)
     */
    ret = LAPACKE_zgeqrf_work( LAPACK_COL_MAJOR, M, rank,
                               u1u2, M, tauU, zbuf, lwork );
    assert( ret == 0 );
    total_flops += FLOPS_ZGEQRF( M, rank );

    /*
     * Concatenate V2 and V1 in v1v2
     *  [ v2^h v2^h v2^h ]
     *  [ 0    v1^h 0    ]
     */
    core_zlrconcatenate_v( transA1, alpha,
                           M1, N1, A,
                               N2, B,
                           offy, v1v2 );

    /*
     * Perform LQ factorization on v1v2 = (L2 Q2)
     */
    ret = LAPACKE_zgelqf_work( LAPACK_COL_MAJOR, rank, N,
                               v1v2, rank, tauV, zbuf, lwork );
    assert(ret == 0);
    total_flops += FLOPS_ZGELQF( rank, N );

    /*
     * Compute R = alpha R1 L2
     */
    u = R + rank * rank;
    v = u + rank * rank;

    memset(R, 0, rank * rank * sizeof(pastix_complex64_t));

    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'U', rank, rank,
                               u1u2, M, R, rank );
    assert(ret == 0);

    cblas_ztrmm(CblasColMajor,
                CblasRight, CblasLower,
                CblasNoTrans, CblasNonUnit,
                rank, rank, CBLAS_SADDR(zone),
                v1v2, rank, R, rank);
    total_flops += FLOPS_ZTRMM( PastixRight, rank, rank );

    if ( lowrank->use_reltol ) {
        /**
         * In relative tolerance, we can choose two solutions:
         *  1) The first one, more conservative, is to compress relatively to
         *  the norm of the final matrix \f$ \alpha A + B \f$. In this kernel, we
         *  exploit the fact that the previous computed product contains all the
         *  information of the final matrix to do it as follow:
         *
         * double norm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', rank, rank,
         *                                    R, rank, NULL );
         * tol = tol * norm;
         *
         *  2) The second solution, less conservative, will allow to reduce the
         *  rank more efficiently. Since A and B have been compressed relatively
         *  to their respective norms, there is no reason to compress the sum
         *  relatively to its own norm, but it is more reasonable to compress it
         *  relatively to the norm of A and B. For example, A-B would be full
         *  with the first criterion, and rank null with the second.
         *  Note that here, we can only have an estimation that once again
         *  reduces the conservation of the criterion.
         *  \f[ || \alpha A + B || <= |\alpha| ||A|| + ||B|| <= |\alpha| ||U_aV_a|| + ||U_bV_b|| \f]
         *
         */
        double normA, normB;
        normA = core_zlrnrm( PastixFrobeniusNorm, transA1,       M1, N1, A );
        normB = core_zlrnrm( PastixFrobeniusNorm, PastixNoTrans, M2, N2, B );
        tol = tol * ( cabs(alpha) * normA + normB );
    }

    /*
     * Compute svd(R) = u sigma v^t
     */
    /* Missing the flops of the u and v generation */
    flops = FLOPS_ZGEQRF( rank, rank ) + FLOPS_ZGELQF( rank, (rank-1) );

    kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_rradd_recompression );
    ret = MYLAPACKE_zgesvd_work( LAPACK_COL_MAJOR, 'S', 'S',
                                 rank, rank, R, rank,
                                 s, u, rank, v, rank,
                                 zbuf, lwork, s + rank );
    if ( ret != 0 ) {
        pastix_print_error( "LAPACKE_zgesvd FAILED" );
    }

    /*
     * Let's compute the new rank of the result
     */
    tmp = v;

    for (i=0; i<rank; i++, tmp+=1){
        double frob_norm;

        /*
         * There are two different stopping criteria for SVD to decide the
         * compression rank:
         *	1) The 2-norm:
         *         Compare the singular values to the threshold
         *      2) The Frobenius norm:
         *         Compare the Frobenius norm of the trailing singular values to
         *         the threshold. Note that we use a reverse accumulation of the
         *         singular values to avoid accuracy issues.
         */
#if defined(PASTIX_SVD_2NORM)
        frob_norm = s[i];
#else
        frob_norm = core_dlassq( rank-i, s + rank - 1, -1 );
#endif

        if (frob_norm < tol) {
            break;
        }
        cblas_zdscal(rank, s[i], tmp, rank);
    }
    new_rank = i;
    kernel_trace_stop_lvl2_rank( flops, new_rank );
    total_flops += flops;

    /*
     * First case: The rank is too big, so we decide to uncompress the result
     */
    if ( new_rank > core_get_rklimit( M, N ) ) {
        pastix_lrblock_t Bbackup = *B;

        core_zlralloc( M, N, -1, B );
        u = B->u;

        /* Uncompress B */
        flops = FLOPS_ZGEMM( M, N, Bbackup.rk );
        kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_uncompress );
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    M, N, Bbackup.rk,
                    CBLAS_SADDR(zone),  Bbackup.u, ldbu,
                                        Bbackup.v, ldbv,
                    CBLAS_SADDR(zzero), u, M );
        kernel_trace_stop_lvl2( flops );
        total_flops += flops;

        /* Add A into it */
        if ( A->rk == -1 ) {
            flops = 2 * M1 * N1;
            kernel_trace_start_lvl2( PastixKernelLvl2_FR_GEMM );
            core_zgeadd( transA1, M1, N1,
                         alpha, A->u, ldau,
                         zone, u + offy * M + offx, M);
            kernel_trace_stop_lvl2( flops );
        }
        else {
            flops = FLOPS_ZGEMM( M1, N1, A->rk );
            kernel_trace_start_lvl2( PastixKernelLvl2_FR_GEMM );
            cblas_zgemm(CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transA1,
                        M1, N1, A->rk,
                        CBLAS_SADDR(alpha), A->u, ldau,
                                            A->v, ldav,
                        CBLAS_SADDR(zone), u + offy * M + offx, M);
            kernel_trace_stop_lvl2( flops );
        }
        total_flops += flops;
        core_zlrfree(&Bbackup);
        memFree_null(zbuf);
        return total_flops;
    }
    else if ( new_rank == 0 ) {
        core_zlrfree(B);
        memFree_null(zbuf);
        return total_flops;
    }

    /*
     * We need to reallocate the buffer to store the new compressed version of B
     * because it wasn't big enough
     */
    ret = core_zlrsze( 0, M, N, B, new_rank, -1, -1 );
    assert( ret != -1 );
    assert( B->rkmax >= new_rank );
    assert( B->rkmax >= B->rk    );

    ldbv = B->rkmax;

    /*
     * Let's now compute the final U = Q1 ([u] sigma)
     *                                     [0]
     */
#if defined(PASTIX_DEBUG_LR_NANCHECK)
    ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, new_rank,
                               0.0, 0.0, B->u, ldbu );
    assert(ret == 0);
    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', rank, new_rank,
                               u, rank, B->u, ldbu );
    assert(ret == 0);
#else
    tmp = B->u;
    for (i=0; i<new_rank; i++, tmp+=ldbu, u+=rank) {
        memcpy(tmp, u,              rank * sizeof(pastix_complex64_t));
        memset(tmp + rank, 0, (M - rank) * sizeof(pastix_complex64_t));
    }
#endif

    flops = FLOPS_ZUNMQR( M, new_rank, minU, PastixLeft )
        +   FLOPS_ZUNMLQ( new_rank, N, minV, PastixRight );
    kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_rradd_computeNewU );
    ret = LAPACKE_zunmqr_work(LAPACK_COL_MAJOR, 'L', 'N',
                              M, new_rank, minU,
                              u1u2, M, tauU,
                              B->u, ldbu,
                              zbuf, lwork);
    assert( ret == 0 );

    /*
     * And the final V^T = [v^t 0 ] Q2
     */
    tmp = B->v;
    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', new_rank, rank,
                               v, rank, B->v, ldbv );
    assert( ret == 0 );

    ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', new_rank, N-rank,
                               0.0, 0.0, tmp + ldbv * rank, ldbv );
    assert( ret == 0 );

    ret = LAPACKE_zunmlq_work(LAPACK_COL_MAJOR, 'R', 'N',
                              new_rank, N, minV,
                              v1v2, rank, tauV,
                              B->v, ldbv,
                              zbuf, lwork);
    assert( ret == 0 );
    kernel_trace_stop_lvl2( flops );
    total_flops += flops;

    (void)ret;
    memFree_null(zbuf);
    return total_flops;
}
