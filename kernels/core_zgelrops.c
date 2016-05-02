/**
 *
 * @file core_zgelrops.c
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

static pastix_complex64_t zzero = 0.;
static pastix_complex64_t zone  = 1.;

//#define WITH_LAPACKE_WORK

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_z_compress_LR - Compresses a dense block into a u v^T LR structure.
 *
 *******************************************************************************
 *
 * @param[in] fL
 *          Pointer to the dense structure of size dimb * dima
 *          Leading dimension is stride
 *
 * @param[out] u
 *          Pointer to the u factor of LR representation of size dimb * rank
 *          Leading dimension is ldu
 *
 * @param[out] v
 *          Pointer to the v factor of LR representation of size dima * rank
 *          Leading dimension is ldv
 *          Note that due to LAPACKE_zgesvd this block is stored transposed
 *
 *
 *******************************************************************************
 *
 * @return
 *          The rank of the compressed structure.
 *
 *******************************************************************************/
int
core_zge2lrx(double tol, pastix_int_t m, pastix_int_t n,
             const pastix_complex64_t *A, pastix_int_t lda,
             pastix_lrblock_t *Alr )
{
    pastix_complex64_t *u, *v, *zwork, *Acpy, ws;
    double             *rwork, *s, tolabs, tolrel;
    pastix_int_t        i, ret, ldu, ldv;
    pastix_int_t        minMN = pastix_imin( m, n );
    pastix_int_t        lwork = -1;
    pastix_int_t        zsize, rsize;

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

    /**
     * Query the workspace needed for the gesvd
     */
    ret = LAPACKE_zgesvd_work( LAPACK_COL_MAJOR, 'S', 'S',
                               m, n, NULL, m,
                               NULL, NULL, ldu, NULL, ldv,
                               &ws, lwork
#if defined(PRECISION_z) || defined(PRECISION_c)
                               , NULL
#endif
                               );

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
#if defined(WITH_LAPACKE_WORK)
    ret = LAPACKE_zlacpy(LAPACK_COL_MAJOR, 'A', m, n,
                         A, lda, Acpy, m );
#else
    ret = LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', m, n,
                              A, lda, Acpy, m );
#endif
    assert( ret == 0 );

#if defined(WITH_LAPACKE_WORK)
    ret = LAPACKE_zgesvd( LAPACK_COL_MAJOR, 'S', 'S',
                          m, n, Acpy, m,
                          s, u, ldu, v, ldv,
                          (double*)zwork );
#else
    ret = LAPACKE_zgesvd_work( LAPACK_COL_MAJOR, 'S', 'S',
                               m, n, Acpy, m,
                               s, u, ldu, v, ldv,
                               zwork, lwork
#if defined(PRECISION_z) || defined(PRECISION_c)
                               , rwork + minMN
#endif
                               );
#endif
    assert(ret == 0);
    if( ret != 0 ){
        errorPrint("SVD Failed\n");
        EXIT(MOD_SOPALIN, INTERNAL_ERR);
    }

    tolrel = tol * s[0];
    tolabs = tol * tol;
    for (i=0; i<minMN; i++, u+=ldu){
        if ( (s[i] >= tolabs) &&
             (s[i] >= tolrel) )
        /* if (s[i] > tol) */
        {
            cblas_zdscal(m, s[i], u, 1);
        }
        else {
            break;
        }
    }
    Alr->rk = i;

    free(zwork);
    return 0;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zge2lr - Convert a full rank matrix in a low rank matrix.
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          The tolerance used as a criterai to eliminate information from the
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
int
core_zge2lr( double tol, pastix_int_t m, pastix_int_t n,
             const pastix_complex64_t *A, pastix_int_t lda,
             pastix_lrblock_t *Alr )
{
    int ret;
    //pastix_complex64_t *tmp;

    Alr->rk = -1;
    Alr->rkmax = pastix_imin( m, n );
    //tmp =  malloc( (m+n) * Alr->rkmax * sizeof(pastix_complex64_t));
    //Alr->u = tmp;
    //Alr->v = tmp + m * Alr->rkmax;
    Alr->u = malloc( m * Alr->rkmax * sizeof(pastix_complex64_t) );
    Alr->v = malloc( n * Alr->rkmax * sizeof(pastix_complex64_t) );

    /**
     * Compress the dense matrix with the temporary space just allocated
     */
    ret = core_zge2lrx( tol, m, n, A, lda, Alr );
    if ( ret != 0 ) {
        free(Alr->u); Alr->u = NULL;
        free(Alr->v); Alr->v = NULL;
        return ret;
    }

    /* The rank should have been set */
    assert(Alr->rk != -1);

    /**
     * It is not interesting to compress, so we store the dense version in Alr
     */
    if ( (Alr->rk * 2) > Alr->rkmax )
    {
        if ( n != Alr->rkmax ) {
            Alr->u = realloc( Alr->u, m * n * sizeof(pastix_complex64_t) );
        }
        free(Alr->v);
        Alr->v = NULL;
        Alr->rk = -1;
        Alr->rkmax = m;
#if defined(WITH_LAPACKE_WORK)
        ret = LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'A', m, n,
                              A, lda, Alr->u, m );
#else
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                                   A, lda, Alr->u, m );
#endif
        assert(ret == 0);
    }
    /**
     * The rank is nul, we free everything
     */
    else if (Alr->rk == 0)
    {
        /**
         * The rank is nul, we free everything
         */
        free(Alr->u);
        free(Alr->v);
        Alr->u = NULL;
        Alr->v = NULL;
        Alr->rkmax = 0;
    }
    /**
     * The rank is non nul, we compress the stored information
     */
    else {
        pastix_complex64_t *u = Alr->u;
        pastix_complex64_t *v = Alr->v;

        Alr->u = malloc( m * Alr->rk * sizeof(pastix_complex64_t) );
        Alr->v = malloc( n * Alr->rk * sizeof(pastix_complex64_t) );

#if defined(WITH_LAPACKE_WORK)
        ret = LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'A', m, Alr->rk,
                                   u, m, Alr->u, m );
        assert(ret == 0);
        ret = LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'A', Alr->rk, n,
                              v, Alr->rkmax, Alr->v, Alr->rk );
#else
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, Alr->rk,
                                   u, m, Alr->u, m );
        assert(ret == 0);
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', Alr->rk, n,
                                   v, Alr->rkmax, Alr->v, Alr->rk );
#endif
        assert(ret == 0);

        free(u); free(v);
        Alr->rkmax = Alr->rk;
    }
    assert( Alr->rk <= Alr->rkmax);

    return 0;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zlr2ge - Convert a low rank matrix into a dense matrix.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          Number of rows of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] n
 *          Number of columns of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] Alr
 *          The low rank matrix to be converted into a dense matrix
 *
 * @param[out] A
 *          The matrix of dimension lda-by-n in which to store the uncompressed
 *          version of Alr.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1, m)
 *
 *******************************************************************************/
int
core_zlr2ge( pastix_int_t m, pastix_int_t n,
             const pastix_lrblock_t *Alr,
             pastix_complex64_t *A, pastix_int_t lda )
{
    int ret;

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
    if ( lda < m ) {
        return -5;
    }
    if ( Alr->rk == -1 ) {
        if (Alr->u == NULL || Alr->v != NULL || Alr->rkmax < m) {
            return -6;
        }
    }
    else if ( Alr->rk != 0){
        if (Alr->u == NULL || Alr->v == NULL) {
            return -6;
        }
    }
#endif

    if ( Alr->rk == -1 ) {
#if defined(WITH_LAPACKE_WORK)
        ret = LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'A', m, n,
                              Alr->u, Alr->rkmax, A, lda );
#else
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                                   Alr->u, Alr->rkmax, A, lda );
#endif
        assert( ret == 0 );
    }
    else if ( Alr->rk == 0 ) {
        ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', m, n,
                                   0., 0., A, lda );
        assert( ret == 0 );
    }
    else {
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    m, n, Alr->rk,
                    CBLAS_SADDR(zone),  Alr->u, m,
                                        Alr->v, Alr->rkmax,
                    CBLAS_SADDR(zzero), A, lda);
    }

    return 0;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_z_add_LR - Adds two LR structure u1 v1^T and (-u2) v2^T into u1 v1^T
 *
 *    u1v1^T + u2v2^T = (u1 u2) (v1 v2)^T
 *    Compute QR decomposition of (u1 u2) = Q1 R1
 *    Compute QR decomposition of (v1 v2) = Q2 R2
 *    Compute SVD of R1 R2^T = u \sigma v^T
 *    Final solution is (Q1 u \sigma^[1/2]) (Q2 v \sigma^[1/2])^T
 *
 *******************************************************************************
 *
 * @param[in, out] u1 v1
 *          LR structure where v1 is stored transposed
 *          u1 factor of size dim_u1 * rank_1 with ld_u1 as leading dimension
 *          v1 factor of size dim_v1 * rank_1 with ld_v1 as leading dimension
 *
 * @param[in] u2 v2
 *          Pointer to the u factor of LR representation of size dimb * rank
 *          Leading dimension is ldu
 *
 * @param[in] x2, y2
 *          Position where u2 v2 is added into u1 v1 (which is larger)
 *
 *
 *******************************************************************************
 *
 * @return
 *          The new rank of u1 v1^T or -1 if ranks are too large for recompression
 *
 *******************************************************************************/
int
core_zrradd( double tol, int transA1, pastix_complex64_t alpha,
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
    double *s, tolabs, tolrel;
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
    ldav = (transA1 == PastixNoTrans) ? A->rkmax : N;
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
             * copying it into B, however the criterai to keep A compressed or
             * not must be based on B dimension, and not on A ones
             */
            u = malloc( M * N * sizeof(pastix_complex64_t));
            B->rk = -1;
            B->rkmax = M;
            B->u = u;
            B->v = NULL;

            if ( M1 != M || N1 != N ) {
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, N,
                                     0., 0., u, M );
            }
            ret = core_zgeadd( PastixNoTrans, M1, N1,
                               alpha, A->u, ldau,
                               0., u + B->rkmax * offy + offx, B->rkmax );
            assert(ret == 0);

            core_zge2lr( tol, M, N, u, M, B );
            free(u);
        }
        else {
            u = malloc( M * A->rk * sizeof(pastix_complex64_t));
            v = malloc( N * A->rk * sizeof(pastix_complex64_t));
            B->rk = A->rk;
            B->rkmax = A->rk;
            B->u = u;
            B->v = v;

            if ( M1 != M ) {
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, B->rk,
                                     0., 0., u, M );
            }
            if ( N1 != N ) {
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', B->rk, N,
                                     0., 0., v, B->rkmax );
            }

#if defined(WITH_LAPACKE_WORK)
            ret = LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'A', M1, A->rk,
                                  A->u, ldau,
                                  u + offx, M );
#else
            ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M1, A->rk,
                                       A->u, ldau,
                                       u + offx, M );
#endif
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
    ret = LAPACKE_zgesvd_work( LAPACK_COL_MAJOR, 'S', 'S',
                               rank, rank, NULL, rank,
                               NULL, NULL, rank, NULL, rank,
                               &querysize, -1
#if defined(PRECISION_z) || defined(PRECISION_c)
                               , NULL
#endif
                               );

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
    //u1u2 = malloc( M * rank * sizeof(pastix_complex64_t));
#if defined(WITH_LAPACKE_WORK)
    LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'A', M, B->rk,
                    B->u, ldbu, u1u2, M );
#else
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, B->rk,
                         B->u, ldbu, u1u2, M );
#endif

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
#if defined(WITH_LAPACKE_WORK)
            ret = LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'A', M1, N1,
                                  A->u, ldau, tmp + offx, M );
#else
            ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M1, N1,
                                       A->u, ldau, tmp + offx, M );
#endif
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
#if defined(WITH_LAPACKE_WORK)
        ret = LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'A', M1, A->rk,
                              A->u, ldau, tmp + offx, M );
#else
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M1, A->rk,
                                   A->u, ldau, tmp + offx, M );
#endif
        assert(ret == 0);
    }

    /**
     * Perform QR factorization on u1u2 = (Q1 R1)
     */
    //tauU = malloc( minU * sizeof(pastix_complex64_t));
#if defined(WITH_LAPACKE_WORK)
    ret = LAPACKE_zgeqrf( LAPACK_COL_MAJOR, M, rank,
                          u1u2, M, tauU );
#else
    ret = LAPACKE_zgeqrf_work( LAPACK_COL_MAJOR, M, rank,
                               u1u2, M, tauU, zbuf, lwork );
#endif
    assert( ret == 0 );

    /**
     * Concatenate V2 and V1 in v1v2
     *  [ v2^h v2^h v2^h ]
     *  [ 0    v1^h 0    ]
     */
    //v1v2 = malloc( N * rank * sizeof(pastix_complex64_t));
#if defined(WITH_LAPACKE_WORK)
    ret = LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'A', B->rk, N,
                          B->v, ldbv, v1v2, rank );
#else
    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', B->rk, N,
                               B->v, ldbv, v1v2, rank );
#endif
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
    //tauV = malloc( minV * sizeof(pastix_complex64_t));
#if defined(WITH_LAPACKE_WORK)
    ret = LAPACKE_zgelqf( LAPACK_COL_MAJOR, rank, N,
                          v1v2, rank, tauV );
#else
    ret = LAPACKE_zgelqf_work( LAPACK_COL_MAJOR, rank, N,
                               v1v2, rank, tauV, zbuf, lwork );
#endif
    assert(ret == 0);
    /**
     * Compute R = alpha R1 L2
     */
    //R = malloc( 3 * rank * rank * sizeof(pastix_complex64_t));
    u = R + rank * rank;
    v = u + rank * rank;

    memset(R, 0, rank * rank * sizeof(pastix_complex64_t));

#if defined(WITH_LAPACKE_WORK)
    ret = LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'U', rank, rank,
                          u1u2, M, R, rank );
#else
    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'U', rank, rank,
                               u1u2, M, R, rank );
#endif
    assert(ret == 0);

    cblas_ztrmm(CblasColMajor,
                CblasRight, CblasLower,
                CblasNoTrans, CblasNonUnit,
                rank, rank, CBLAS_SADDR(zone),
                v1v2, rank, R, rank);

    /**
     * Compute svd(R) = u \sigma v^t
     */
#if defined(WITH_LAPACKE_WORK)
    ret = LAPACKE_zgesvd( CblasColMajor, 'S', 'S',
                          rank, rank, R, rank,
                          s, u, rank, v, rank,
                          (double*)zbuf);
#else
    ret = LAPACKE_zgesvd_work( CblasColMajor, 'S', 'S',
                               rank, rank, R, rank,
                               s, u, rank, v, rank,
                               zbuf, lwork
#if defined(PRECISION_z) || defined(PRECISION_c)
                               , s + rank
#endif
                               );
#endif
    assert(ret == 0);
    if (ret != 0) {
        errorPrint("LAPACKE_zgesvd FAILED");
        EXIT(MOD_SOPALIN, INTERNAL_ERR);
    }

    /**
     * Let's compute the new rank of the result
     */
    tmp = u;
    tolrel = tol * s[0];
    tolabs = tol * tol;
    for (i=0; i<rank; i++, tmp+=rank){
        if ( (s[i] >= tolabs) &&
             (s[i] >= tolrel) )
        /* if (s[i] > tol) */
        {
            cblas_zdscal(rank, s[i], tmp, 1);
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
        fprintf(stderr, "The new rank became too large, we uncompress!!!!\n");
        zbuf = realloc( zbuf, M * N * sizeof(pastix_complex64_t) );

        /* Uncompress B */
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    M, N, B->rk,
                    CBLAS_SADDR(zone),  B->u, ldbu,
                                        B->v, ldbv,
                    CBLAS_SADDR(zzero), zbuf, M);

        /* Add A into it */
        if ( A->rk == -1 ) {
            core_zgeadd( transA1, M1, N1,
                         alpha, A->u, ldau,
                         zone, zbuf + offy * M + offx, M);
        }
        else {
            cblas_zgemm(CblasColMajor, CblasNoTrans, transA1,
                        M1, N1, A->rk,
                        CBLAS_SADDR(alpha), A->u, ldau,
                                            A->v, ldav,
                        CBLAS_SADDR(zone), zbuf + offy * M + offx, M);
        }
        free(B->u);
        free(B->v);
        B->rk = -1;
        B->rkmax = M;
        B->u = zbuf;
        B->v = NULL;
        assert( B->rk <= B->rkmax);
        return 0;
    }
    else if ( new_rank == 0 ) {
        free(B->u);
        free(B->v);
        B->rk = 0;
        B->rkmax = 0;
        B->u = NULL;
        B->v = NULL;
        free(zbuf);
        return 0;
    }

    B->rk = new_rank;
    assert(rank <= M && rank <= N);

    /**
     * We need to reallocate the buffer to store the new compressed version of B
     * because it wasn't big enough
     */
    if ( new_rank > B->rkmax ) {
        /**
         * We use a temporary buffer to allow pointer arithmetic
         */
        /* pastix_complex64_t *uv = malloc( new_rank * (M+N) * sizeof(pastix_complex64_t)); */

        /* free(B->u); */
        /* B->u = uv; */
        /* B->v = uv + M * new_rank; */
        /* B->rkmax = new_rank; */

        free(B->u);
        free(B->v);
        B->u = malloc( new_rank * M * sizeof(pastix_complex64_t));
        B->v = malloc( new_rank * N * sizeof(pastix_complex64_t));
        ldbv = new_rank;
#if defined(WITH_LAPACKE_WORK)
        LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, new_rank,
                             NAN, NAN, B->u, ldbu );
        LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', new_rank, N,
                             NAN, NAN, B->v, ldbv );
#endif

        B->rkmax = new_rank;
    }
    assert( B->rkmax >= new_rank );
    assert( B->rkmax >= B->rk    );

    /**
     * Let's now compute the final U = Q1 ([u] \sigma)
     *                                     [0]
     */
#if defined(WITH_LAPACKE_WORK)
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

#if defined(WITH_LAPACKE_WORK)
    ret = LAPACKE_zunmqr(LAPACK_COL_MAJOR, 'L', 'N',
                              M, new_rank, minU,
                              u1u2, M, tauU,
                              B->u, ldbu);
#else
    ret = LAPACKE_zunmqr_work(LAPACK_COL_MAJOR, 'L', 'N',
                              M, new_rank, minU,
                              u1u2, M, tauU,
                              B->u, ldbu,
                              zbuf, lwork);
#endif
    assert( ret == 0 );

    /**
     * And the final V^T = [v^t 0 ] Q2
     */
    tmp = B->v;
#if defined(WITH_LAPACKE_WORK)
    ret = LAPACKE_zlacpy( LAPACK_COL_MAJOR, 'A', new_rank, rank,
                         v, rank, B->v, ldbv );
#else
    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', new_rank, rank,
                               v, rank, B->v, ldbv );
#endif
    assert( ret == 0 );

    ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', new_rank, N-rank,
                               0., 0., tmp + ldbv * rank, ldbv );
    assert( ret == 0 );

#if defined(WITH_LAPACKE_WORK)
    ret = LAPACKE_zunmlq(LAPACK_COL_MAJOR, 'R', 'N',
                         new_rank, N, minV,
                         v1v2, rank, tauV,
                         B->v, ldbv);
#else
    ret = LAPACKE_zunmlq_work(LAPACK_COL_MAJOR, 'R', 'N',
                              new_rank, N, minV,
                              v1v2, rank, tauV,
                              B->v, ldbv,
                              zbuf, lwork);
#endif
    assert( ret == 0 );

    free(zbuf);
    //free(R);
    return new_rank;
}

int
core_zgradd( double tol, pastix_complex64_t alpha,
             pastix_int_t M1, pastix_int_t N1, pastix_complex64_t *A, pastix_int_t lda,
             pastix_int_t M2, pastix_int_t N2, pastix_lrblock_t *B,
             pastix_int_t offx, pastix_int_t offy)
{
    pastix_lrblock_t lrA;
    pastix_int_t rmax = pastix_imin( M2, N2 );
    pastix_int_t rank, ldub;

    assert( B->rk <= B->rkmax);

    if ( B->rk == -1 ) {
        pastix_complex64_t *tmp = B->u;
        ldub = B->rkmax;
        tmp += ldub * offy + offx;
        core_zgeadd( CblasNoTrans, M1, N1,
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
                     1.,    work + M2 * offy + offx, M2 );

        core_zge2lr( tol, M2, N2, work, M2, B );
        rank = B->rk;
        free(work);
    }
    /**
     * We consider the A matrix as Id * A or A *Id
     */
    else {
        lrA.rk = -1;
        lrA.rkmax = lda;
        lrA.u = A;
        lrA.v = NULL;
        rank = core_zrradd( tol, PastixNoTrans, alpha,
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
 * @ingroup pastix_kernel
 *
 * core_zlrm2 - Computes the product of two low rank matrices and returns the result in AB
 *
 *******************************************************************************/
int core_zlrm2( int transA, int transB,
                pastix_int_t M, pastix_int_t N, pastix_int_t K,
                const pastix_lrblock_t *A,
                const pastix_lrblock_t *B,
                pastix_lrblock_t *AB,
                pastix_complex64_t *work,
                pastix_int_t ldwork )
{
    pastix_int_t ldau, ldav, ldbu, ldbv;
    int transV = PastixNoTrans;

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
            /**
             * Let's compute A * B' = Au Av^h (Bu Bv^h)' with the smallest ws
             */
            if ( (A->rk * N) <= (B->rk * M) ) {
                /**
                 *    ABu = Au
                 *    ABv = (Av^h Bv^h') * Bu'
                 */
                assert( (A->rk * ( N + B->rk )) <= ldwork );
                AB->rk = A->rk;
                AB->rkmax = A->rk;
                AB->u = A->u;
                AB->v = work + A->rk * B->rk;

                cblas_zgemm( CblasColMajor, CblasNoTrans, transB,
                             A->rk, B->rk, K,
                             CBLAS_SADDR(zone),  A->v, ldav,
                                                 B->v, ldbv,
                             CBLAS_SADDR(zzero), work, A->rk );

                cblas_zgemm( CblasColMajor, CblasNoTrans, transB,
                             A->rk, N, B->rk,
                             CBLAS_SADDR(zone),  work,  A->rk,
                                                 B->u,  ldbu,
                             CBLAS_SADDR(zzero), AB->v, AB->rkmax );
            }
            else {
                /**
                 *    ABu = Au * (Av^h Bv^h')
                 *    ABv = Bu'
                 */
                assert( (B->rk * ( M + A->rk )) <= ldwork );
                AB->rk = B->rk;
                AB->rkmax = B->rk;
                AB->u = work + A->rk * B->rk;
                AB->v = B->u;

                cblas_zgemm( CblasColMajor, CblasNoTrans, transB,
                             A->rk, B->rk, K,
                             CBLAS_SADDR(zone),  A->v, ldav,
                                                 B->v, ldbv,
                             CBLAS_SADDR(zzero), work, A->rk );

                cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                             M, B->rk, A->rk,
                             CBLAS_SADDR(zone),  A->u,  ldau,
                                                 work,  A->rk,
                             CBLAS_SADDR(zzero), AB->u, M );

                transV = transB;
            }
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

            cblas_zgemm( CblasColMajor, CblasNoTrans, transB,
                         A->rk, N, K,
                         CBLAS_SADDR(zone),  A->v,  ldav,
                                             B->u,  ldbu,
                         CBLAS_SADDR(zzero), AB->v, AB->rkmax );
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

            cblas_zgemm( CblasColMajor, CblasNoTrans, transB,
                         M, B->rk, K,
                         CBLAS_SADDR(zone),  A->u,  ldau,
                                             B->v,  ldbv,
                         CBLAS_SADDR(zzero), AB->u, M );

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
            if ( 2*K < pastix_imin( M, N ) ) {
                AB->rk = K;
                AB->rkmax = B->rkmax;
                AB->u = A->u;
                AB->v = B->u;
                transV = transB;
            }
            else {
                assert( (M * N) <= ldwork );
                AB->rk = -1;
                AB->rkmax = M;
                AB->u = work;
                AB->v = NULL;

                cblas_zgemm( CblasColMajor, CblasNoTrans, transB,
                             M, N, K,
                             CBLAS_SADDR(zone),  A->u, ldau,
                                                 B->u, ldbu,
                             CBLAS_SADDR(zzero), work, M );
            }
        }
    }
    assert( AB->rk <= AB->rkmax);
    return transV;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zlrmm - A * B + C with three Low rank matrices
 *
 *******************************************************************************/
int
core_zlrmm( double tol, int transA, int transB,
            pastix_int_t M, pastix_int_t N, pastix_int_t K,
            pastix_int_t Cm, pastix_int_t Cn,
            pastix_int_t offx, pastix_int_t offy,
            pastix_complex64_t alpha, const pastix_lrblock_t *A,
                                      const pastix_lrblock_t *B,
            pastix_complex64_t beta,        pastix_lrblock_t *C,
            pastix_complex64_t *work, pastix_int_t ldwork )
{
    pastix_complex64_t *tmp = NULL;
    pastix_lrblock_t AB;
    pastix_int_t ldabu, ldabv, ldcu, ldcv;
    pastix_int_t required = 0;
    int transV;
    int allocated = 0;

    assert(transA == PastixNoTrans);
    assert(transB != PastixNoTrans);
    assert( A->rk <= A->rkmax);
    assert( B->rk <= B->rkmax);
    assert( C->rk <= C->rkmax);

    /* Quick return if multiplication by 0 */
    if ( A->rk == 0 || B->rk == 0 ) {
        return 0;
    }

    if ( A->rk != -1 ) {
        if ( B->rk != -1 ) {
            required = pastix_imin( A->rk * ( N + B->rk ),
                                    B->rk * ( M + A->rk ) );
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

    transV = core_zlrm2( transA, transB, M, N, K,
                         A, B, &AB, tmp, required );

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
            cblas_zgemm( CblasColMajor, CblasNoTrans, transV,
                         M, N, AB.rk,
                         CBLAS_SADDR(alpha), AB.u, ldabu,
                                             AB.v, ldabv,
                         CBLAS_SADDR(beta),  Cptr, ldcu );
        }
    }
     /**
     * The destination matrix is low rank
     */
    else {
        if ( AB.rk == -1 ) {
            assert(beta == 1.);
            core_zgradd( tol, alpha,
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
                    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                                 Cm, Cn, C->rk,
                                 CBLAS_SADDR(beta),  C->u, ldcu,
                                 C->v, ldcv,
                                 CBLAS_SADDR(zzero), work, Cm );
                }
                else{
                    memset(work, 0, Cm * Cn * sizeof(pastix_complex64_t) );
                }

                /* Add A*B */
                cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                             M, N, AB.rk,
                             CBLAS_SADDR(alpha), AB.u, M,
                                                 AB.v, AB.rkmax,
                             CBLAS_SADDR(zone), work + Cm * offy + offx, Cm );

                core_zge2lr( tol, Cm, Cn, work, Cm, C );
                free(work);
            }
            else {
                /* Need to handle correctly this case */
                core_zrradd( tol, transV, alpha,
                             M, N, &AB,
                             Cm, Cn, C,
                             offx, offy );
            }
        }
    }

    if ( allocated ) {
        free(tmp);
    }

    assert( C->rk <= C->rkmax);
    return 0;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zlrmge - A * B + C with A, and B Low rank matrices, and C full rank
 *
 *******************************************************************************/
int
core_zlrmge( double tol, int transA, int transB,
             pastix_int_t M, pastix_int_t N, pastix_int_t K,
             pastix_complex64_t alpha, const pastix_lrblock_t *A,
                                       const pastix_lrblock_t *B,
             pastix_complex64_t beta, pastix_complex64_t *C, int ldc,
             pastix_complex64_t *work, pastix_int_t ldwork )
{
    pastix_lrblock_t lrC;

    lrC.rk = -1;
    lrC.rkmax = ldc;
    lrC.u = C;
    lrC.v = NULL;

    core_zlrmm( tol, transA, transB, M, N, K,
                M, N, 0, 0,
                alpha, A, B, beta, &lrC,
                work, ldwork );

    return 0;
}
