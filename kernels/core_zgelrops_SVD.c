/**
 *
 * @file core_zgelrops_SVD.c
 *
 *  PaStiX kernel routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Gregoire Pichon
 * @date 2016-23-03
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <cblas.h>
#include <lapacke.h>
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "z_nan_check.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static pastix_complex64_t zone  =  1.;
static pastix_complex64_t zzero =  0.;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zge2lrx - Computes a SVD with truncation.
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          The tolerance used as a criteria to eliminate information from the
 *          full rank matrix
 *
 * @param[in] m
 *          Number of rows of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] n
 *          Number of columns of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] A
 *          The matrix of dimension lda-by-n that need to be compressed
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
 * @return
 *          This routine will return the rank of A
 *
 *******************************************************************************/
int
core_zge2lrx(double tol, pastix_int_t m, pastix_int_t n,
             const pastix_complex64_t *A, pastix_int_t lda,
             pastix_lrblock_t *Alr )
{
    pastix_complex64_t *u, *v, *zwork, *Acpy, ws;
    double             *rwork, *s;
    /* double              tolabs, tolrel; */
    pastix_int_t        i, ret, ldu, ldv;
    pastix_int_t        minMN = pastix_imin( m, n );
    pastix_int_t        lwork = -1;
    pastix_int_t        zsize, rsize;
    double              norm, relative_tolerance;

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
    if ( (Alr->u == NULL) || (Alr->v == NULL) || (Alr->rkmax < minMN) ) {
        return -6;
    }
#endif

    u = Alr->u;
    v = Alr->v;
    ldu = m;
    ldv = Alr->rkmax;

    norm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                A, lda, NULL );
    relative_tolerance = tol * norm;

    /**
     * Query the workspace needed for the gesvd
     */
#if defined(PASTIX_DEBUG_LR_NANCHECK)
    ws = minMN;
#else
    ret = MYLAPACKE_zgesvd_work( LAPACK_COL_MAJOR, 'S', 'S',
                                 m, n, NULL, m,
                                 NULL, NULL, ldu, NULL, ldv,
                                 &ws, lwork, NULL );
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

    /**
     * Backup the original matrix before to overwrite it with the SVD
     */
    ret = LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', m, n,
                              A, lda, Acpy, m );
    assert( ret == 0 );

    ret = MYLAPACKE_zgesvd_work( LAPACK_COL_MAJOR, 'S', 'S',
                                 m, n, Acpy, m,
                                 s, u, ldu, v, ldv,
                                 zwork, lwork, rwork + minMN );
    assert(ret == 0);
    if( ret != 0 ){
        errorPrint("SVD Failed\n");
        EXIT(MOD_SOPALIN, INTERNAL_ERR);
    }

    for (i=0; i<minMN; i++, v+=1){
        if (s[i] > relative_tolerance)
        {
            cblas_zdscal(n, s[i], v, ldv);
        }
        else {
            break;
        }
    }
    Alr->rk = i;

    memFree_null(zwork);
    return i;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zge2lr_SVD - Convert a full rank matrix in a low rank matrix, using SVD.
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          The tolerance used as a criteria to eliminate information from the
 *          full rank matrix
 *
 * @param[in] m
 *          Number of rows of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] n
 *          Number of columns of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] A
 *          The matrix of dimension lda-by-n that need to be compressed
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1, m)
 *
 * @param[out] Alr
 *          The low rank matrix structure that will store the low rank
 *          representation of A
 *
 *******************************************************************************/
void
core_zge2lr_SVD( double tol, pastix_int_t m, pastix_int_t n,
                 const pastix_complex64_t *A, pastix_int_t lda,
                 pastix_lrblock_t *Alr )
{
    int ret;
    /**
     * Allocate a temorary Low rank matrix
     */
    core_zlralloc( m, n, pastix_imin( m, n ), Alr );

    /**
     * Compress the dense matrix with the temporary space just allocated
     */
    ret = core_zge2lrx( tol, m, n, A, lda, Alr );

    if ( ret < 0 ) {
        core_zlrfree( Alr );
    }

    /**
     * Resize the space used by the low rank matrix
     */
    ret = core_zlrsze( 1, m, n, Alr, ret, -1 );

    /**
     * It was not interesting to compress, so we store the dense version in Alr
     */
    if (ret == -1) {
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                                   A, lda, Alr->u, Alr->rkmax );
        assert(ret == 0);
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zrradd_SVD - Adds two LR structures A=(-u1) v1^T and B=u2 v2^T into u2 v2^T
 *
 *    u2v2^T - u1v1^T = (u2 u1) (v2 v1)^T
 *    Compute QR decomposition of (u2 u1) = Q1 R1
 *    Compute QR decomposition of (v2 v1) = Q2 R2
 *    Compute SVD of R1 R2^T = u sigma v^T
 *    Final solution is (Q1 u sigma^[1/2]) (Q2 v sigma^[1/2])^T
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          The absolute tolerance criteria
 *
 * @param[in] transA1
 *         @arg CblasNoTrans   :  No transpose, op( A ) = A;
 *         @arg CblasTrans     :  Transpose, op( A ) = A';
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
 * @return
 *          The new rank of u2 v2^T or -1 if ranks are too large for recompression
 *
 *******************************************************************************/
int
core_zrradd_SVD( double tol, pastix_trans_t transA1, pastix_complex64_t alpha,
                 pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                 pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                 pastix_int_t offx, pastix_int_t offy)
{
    pastix_int_t rank, M, N, minU, minV;
    pastix_int_t i, ret, lwork, new_rank;
    pastix_int_t ldau, ldav, ldbu, ldbv;
    pastix_complex64_t *u1u2, *v1v2, *R, *u, *v;
    pastix_complex64_t *tmp, *zbuf, *tauU, *tauV;
    pastix_complex64_t  querysize;
    double *s;
    double norm, relative_tolerance;
    size_t wzsize, wdsize;

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
        errorPrint("Dimensions are not correct");
        assert(0 /* Incorrect dimensions */);
        return -1;
    }

    /**
     * A is rank null, nothing to do
     */
    if (A->rk == 0) {
        return rank;
    }

    ldau = (A->rk == -1) ? A->rkmax : M1;
    ldav = (transA1 == PastixNoTrans) ? A->rkmax : N1;
    ldbu = M;
    ldbv = B->rkmax;

    /**
     * Let's handle case where B is a null matrix
     *   B = alpha A
     */
    if (B->rk == 0) {
        if ( A->rk == -1 ) {
            /**
             * TODO: This case can be improved by compressing A, and then
             * copying it into B, however the criteria to keep A compressed or
             * not must be based on B dimension, and not on A ones
             */
            MALLOC_INTERN( u, M * N, pastix_complex64_t );

            if ( M1 != M || N1 != N ) {
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, N,
                                     0., 0., u, M );
            }
            ret = core_zgeadd( PastixNoTrans, M1, N1,
                               alpha, A->u, ldau,
                               0., u + M * offy + offx, M );
            assert(ret == 0);

            core_zge2lr_SVD( tol, M, N, u, M, B );
            memFree_null(u);
        }
        else {
            core_zlralloc( M, N, A->rkmax, B );
            u = B->u;
            v = B->v;
            B->rk = A->rk;

            if ( M1 != M ) {
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, B->rk,
                                     0., 0., u, M );
            }
            if ( N1 != N ) {
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', B->rk, N,
                                     0., 0., v, B->rkmax );
            }

            ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M1, A->rk,
                                       A->u, ldau,
                                       u + offx, M );
            assert(ret == 0);

            ret = core_zgeadd( transA1, A->rk, N1,
                               alpha, A->v, ldav,
                               0., v + B->rkmax * offy, B->rkmax );
            assert(ret == 0);
        }
        assert( B->rk <= B->rkmax);
        return 0;
    }

    /**
     * The rank is too big, let's try to compress
     */
    if ( rank > pastix_imin( M, N ) ) {
        assert(0);
    }

    /**
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
    ret = MYLAPACKE_zgesvd_work( LAPACK_COL_MAJOR, 'S', 'S',
                                 rank, rank, NULL, rank,
                                 NULL, NULL, rank, NULL, rank,
                                 &querysize, -1, NULL);
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

    /**
     * Concatenate U2 and U1 in u1u2
     *  [ u2  0  ]
     *  [ u2  u1 ]
     *  [ u2  0  ]
     */
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, B->rk,
                         B->u, ldbu, u1u2, M );

    tmp = u1u2 + B->rk * M;
    if ( A->rk == -1 ) {
        /**
         * A is full of rank M1, so A will be integrated into v1v2
         */
        if ( M1 < N1 ) {
            if (M1 != M2) {
                /* Set to 0 */
                memset(tmp, 0, M * M1 * sizeof(pastix_complex64_t));

                /* Set diagonal */
                tmp += offx;
                for (i=0; i<M1; i++, tmp += M+1) {
                    *tmp = 1.;
                }
            }
            else {
                assert( offx == 0 );
                ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, M1,
                                           0., 1., tmp, M );
                assert( ret == 0 );
            }
        }
        else {
            /**
             * A is full of rank N1, so A is integrated into u1u2
             */
            if (M1 != M) {
                memset(tmp, 0, M * N1 * sizeof(pastix_complex64_t));
            }
            ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M1, N1,
                                       A->u, ldau, tmp + offx, M );
            assert(ret == 0);
        }
    }
    /**
     * A is low rank of rank A->rk
     */
    else {
        if (M1 != M) {
            memset(tmp, 0, M * A->rk * sizeof(pastix_complex64_t));
        }
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M1, A->rk,
                                   A->u, ldau, tmp + offx, M );
        assert(ret == 0);
    }

    /**
     * Perform QR factorization on u1u2 = (Q1 R1)
     */
    ret = LAPACKE_zgeqrf_work( LAPACK_COL_MAJOR, M, rank,
                               u1u2, M, tauU, zbuf, lwork );
    assert( ret == 0 );

    /**
     * Concatenate V2 and V1 in v1v2
     *  [ v2^h v2^h v2^h ]
     *  [ 0    v1^h 0    ]
     */
    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', B->rk, N,
                               B->v, ldbv, v1v2, rank );
    assert(ret == 0);

    tmp = v1v2 + B->rk;
    if ( A->rk == -1 ) {
        assert( transA1 == PastixNoTrans );
        /**
         * A is full of rank M1, so it is integrated into v1v2
         */
        if ( M1 < N1 ) {
            if (N1 != N) {
                ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M1, N,
                                           0., 0., tmp, rank );
                assert( ret == 0 );
            }
            core_zgeadd( PastixNoTrans, M1, N1,
                         alpha, A->u, ldau,
                         0., tmp + offy * rank, rank );
        }
        /**
         * A is full of rank N1, so it has been integrated into u1u2
         */
        else {
            if (N1 != N2) {
                /* Set to 0 */
                ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', N1, N,
                                           0., 0., tmp, rank );
                assert(ret == 0);

                /* Set diagonal */
                tmp += offy * rank;
                for (i=0; i<N1; i++, tmp += rank+1) {
                    *tmp = alpha;
                }
            }
            else {
                assert( offy == 0 );
                ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', N1, N,
                                           0., alpha, tmp + offy * rank, rank );
                assert( ret == 0 );
            }
        }
    }
    /**
     * A is low rank of rank A->rk
     */
    else {
        if (N1 != N) {
            ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', A->rk, N,
                                       0., 0., tmp, rank );
            assert(ret == 0);
        }
        core_zgeadd( transA1, A->rk, N1,
                     alpha, A->v,              ldav,
                        0., tmp + offy * rank, rank );
    }

    /**
     * Perform LQ factorization on v1v2 = (L2 Q2)
     */
    ret = LAPACKE_zgelqf_work( LAPACK_COL_MAJOR, rank, N,
                               v1v2, rank, tauV, zbuf, lwork );
    assert(ret == 0);
    /**
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

    norm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', rank, rank,
                               R, rank, NULL );
    relative_tolerance = tol * norm;

    /**
     * Compute svd(R) = u sigma v^t
     */
    ret = MYLAPACKE_zgesvd_work( CblasColMajor, 'S', 'S',
                                 rank, rank, R, rank,
                                 s, u, rank, v, rank,
                                 zbuf, lwork, s + rank );
    assert(ret == 0);
    if (ret != 0) {
        errorPrint("LAPACKE_zgesvd FAILED");
        EXIT(MOD_SOPALIN, INTERNAL_ERR);
    }

    /**
     * Let's compute the new rank of the result
     */
    tmp = v;

    for (i=0; i<rank; i++, tmp+=1){
        if (s[i] > relative_tolerance)
        {
            cblas_zdscal(rank, s[i], tmp, rank);
        }
        else {
            break;
        }
    }
    new_rank = i;

    /**
     * First case: The rank is too big, so we decide to uncompress the result
     */
    if ( new_rank*2 > pastix_imin( M, N ) ) {
        pastix_lrblock_t Bbackup = *B;

        core_zlralloc( M, N, -1, B );
        u = B->u;

        /* Uncompress B */
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    M, N, Bbackup.rk,
                    CBLAS_SADDR(zone),  Bbackup.u, ldbu,
                                        Bbackup.v, ldbv,
                    CBLAS_SADDR(zzero), u, M );

        /* Add A into it */
        if ( A->rk == -1 ) {
            core_zgeadd( transA1, M1, N1,
                         alpha, A->u, ldau,
                         zone, u + offy * M + offx, M);
        }
        else {
            cblas_zgemm(CblasColMajor, CblasNoTrans, transA1,
                        M1, N1, A->rk,
                        CBLAS_SADDR(alpha), A->u, ldau,
                                            A->v, ldav,
                        CBLAS_SADDR(zone), u + offy * M + offx, M);
        }
        core_zlrfree(&Bbackup);
        memFree_null(zbuf);
        return 0;
    }
    else if ( new_rank == 0 ) {
        core_zlrfree(B);
        memFree_null(zbuf);
        return 0;
    }

    /**
     * We need to reallocate the buffer to store the new compressed version of B
     * because it wasn't big enough
     */
    ret = core_zlrsze( 0, M, N, B, new_rank, -1 );
    assert( ret != -1 );
    assert( B->rkmax >= new_rank );
    assert( B->rkmax >= B->rk    );

    ldbv = B->rkmax;

    /**
     * Let's now compute the final U = Q1 ([u] sigma)
     *                                     [0]
     */
#if defined(PASTIX_DEBUG_LR_NANCHECK)
    ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, new_rank,
                               0., 0., B->u, ldbu );
    assert(ret == 0);
    ret = LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'A', rank, new_rank,
                          u, rank, B->u, ldbu );
    assert(ret == 0);
#else
    tmp = B->u;
    for (i=0; i<new_rank; i++, tmp+=ldbu, u+=rank) {
        memcpy(tmp, u,              rank * sizeof(pastix_complex64_t));
        memset(tmp + rank, 0, (M - rank) * sizeof(pastix_complex64_t));
    }
#endif

    ret = LAPACKE_zunmqr_work(LAPACK_COL_MAJOR, 'L', 'N',
                              M, new_rank, minU,
                              u1u2, M, tauU,
                              B->u, ldbu,
                              zbuf, lwork);
    assert( ret == 0 );

    /**
     * And the final V^T = [v^t 0 ] Q2
     */
    tmp = B->v;
    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', new_rank, rank,
                               v, rank, B->v, ldbv );
    assert( ret == 0 );

    ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', new_rank, N-rank,
                               0., 0., tmp + ldbv * rank, ldbv );
    assert( ret == 0 );

    ret = LAPACKE_zunmlq_work(LAPACK_COL_MAJOR, 'R', 'N',
                              new_rank, N, minV,
                              v1v2, rank, tauV,
                              B->v, ldbv,
                              zbuf, lwork);
    assert( ret == 0 );

    memFree_null(zbuf);
    return new_rank;
}

/* Interfaces to transform pastix_complex64_t into void */
/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_kernel
 *
 * core_zge2lr_SVD_interface - Interface to core_zge2lr_SVD
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          The tolerance used as a criteria to eliminate information from the
 *          full rank matrix
 *
 * @param[in] m
 *          Number of rows of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] n
 *          Number of columns of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] Aptr
 *          The matrix of dimension lda-by-n that need to be compressed
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1, m)
 *
 * @param[out] Alr
 *          The low rank matrix structure that will store the low rank
 *          representation of A
 *
 *******************************************************************************/
void core_zge2lr_SVD_interface( pastix_fixdbl_t tol, pastix_int_t m, pastix_int_t n,
                                const void *Aptr, pastix_int_t lda,
                                void *Alr )
{
    const pastix_complex64_t *A = (const pastix_complex64_t *) Aptr;
    core_zge2lr_SVD( tol, m, n, A, lda, Alr );
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_devkernel
 *
 * core_zrradd_SVD_interface - Interface to core_zrradd_SVD
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          The absolute tolerance criteria
 *
 * @param[in] transA1
 *         @arg CblasNoTrans   :  No transpose, op( A ) = A;
 *         @arg CblasTrans     :  Transpose, op( A ) = A';
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
 * @return
 *          The new rank of u2 v2^T or -1 if ranks are too large for recompression
 *
 *******************************************************************************/
int core_zrradd_SVD_interface( pastix_fixdbl_t tol, pastix_trans_t transA1, const void *alphaptr,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy)
{
    const pastix_complex64_t *alpha = (const pastix_complex64_t *) alphaptr;
    return core_zrradd_SVD( tol, transA1, *alpha,
                            M1, N1, A,
                            M2, N2, B,
                            offx, offy );
}
