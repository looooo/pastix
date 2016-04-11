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
static pastix_complex64_t mzone = -1.;

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
pastix_int_t
core_z_compress_LR(double tol, pastix_int_t m, pastix_int_t n,
                   const pastix_complex64_t *A, pastix_int_t lda,
                   pastix_complex64_t *u, pastix_int_t ldu,
                   pastix_complex64_t *v, pastix_int_t ldv)
{
    pastix_complex64_t *zwork, *Acpy, ws;
    double             *rwork, *s;
    pastix_int_t        ret;
    pastix_int_t        i;
    pastix_int_t        minMN = pastix_imin( m, n );
    pastix_int_t        rank = minMN;
    pastix_int_t        lwork = -1;
    pastix_int_t        zsize, rsize;

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

    /* TODO: remove this copy, we can erase the matrix in most cases */
    Acpy = zwork + lwork;
    s    = rwork;

    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', m, n,
                        A, lda, Acpy, m );

    ret = LAPACKE_zgesvd_work( LAPACK_COL_MAJOR, 'S', 'S',
                               m, n, Acpy, m,
                               s, u, ldu, v, ldv,
                               zwork, lwork
#if defined(PRECISION_z) || defined(PRECISION_c)
                               , rwork + minMN
#endif
                               );

    assert(ret == 0);
    if( ret != 0 ){
        errorPrint("SVD Failed\n");
        EXIT(MOD_SOPALIN, INTERNAL_ERR);
    }

    for (i=0; i<minMN; i++){
        if (s[i] < tol) { /// s[0] < tolerance || s[i] < 0.5*tolerance){
            rank = i+1;
            break;
        }
    }

    /* Scale u with the singular values */
    for (i=0; i<rank; i++, u += ldu){
        cblas_zscal(m, CBLAS_SADDR(s[i]), u, 1);
    }

    free(zwork);
    return rank;
}

void core_zge2lr( double tol, pastix_int_t m, pastix_int_t n,
                  const pastix_complex64_t *A, pastix_int_t lda,
                  pastix_lrblock_t *Alr )
{
    Alr->rkmax = pastix_imin( m, n );
    Alr->u = malloc( m * Alr->rkmax * sizeof(pastix_complex64_t));
    Alr->v = malloc( n * Alr->rkmax * sizeof(pastix_complex64_t));

    Alr->rk = core_z_compress_LR( tol, m, n, A, lda,
                                  Alr->u, m, Alr->v, Alr->rkmax );

    /* Could not compress */
    if ( Alr->rk == -1 ) {
        if ( n != Alr->rkmax ) {
            Alr->u = realloc( Alr->u, m * n * sizeof(pastix_complex64_t) );
        }
        free( Alr->v );
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                             A, lda, Alr->u, m );
    }
}

void core_zlr2ge( pastix_int_t m, pastix_int_t n,
                  pastix_lrblock_t *Alr,
                  pastix_complex64_t *A, pastix_int_t lda )
{
    if ( Alr->rk == -1 ) {
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                             Alr->u, m, A, lda );
        free(Alr->u);
    }
    else {
        core_z_uncompress_LR( m, n, Alr->rk,
                              Alr->u, m, Alr->v, Alr->rkmax,
                              A, lda );
        free(Alr->u);
        free(Alr->v);
        Alr->rk = -1;
        Alr->rkmax = -1;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_z_uncompress_LR - Uncompresses a u v^T LR structure into a dense block.
 *
 *******************************************************************************
 *
 * @param[out] fL
 *          Pointer to the dense structure of size dimb * dima
 *          Leading dimension is stride
 *
 * @param[in] u
 *          Pointer to the u factor of LR representation of size dimb * rank
 *          Leading dimension is ldu
 *
 * @param[in] v
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
void
core_z_uncompress_LR( pastix_int_t m, pastix_int_t n, pastix_int_t rank,
                      const pastix_complex64_t *u, pastix_int_t ldu,
                      const pastix_complex64_t *v, pastix_int_t ldv,
                            pastix_complex64_t *A, pastix_int_t lda )
{
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                m, n, rank,
                CBLAS_SADDR(zone),  u, ldu,
                                    v, ldv,
                CBLAS_SADDR(zzero), A, lda);
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
pastix_int_t
core_z_add_LR2(double tol, pastix_complex64_t alpha,
               pastix_int_t M1, pastix_int_t N1, pastix_int_t r1,
               const pastix_complex64_t *u1, pastix_int_t ldu1,
               const pastix_complex64_t *v1, pastix_int_t ldv1,
               pastix_int_t M2, pastix_int_t N2, pastix_int_t r2,
               pastix_complex64_t *u2, pastix_int_t ldu2,
               pastix_complex64_t *v2, pastix_int_t ldv2,
               pastix_int_t offx, pastix_int_t offy)
{
    pastix_int_t rank  = r1 + r2;
    pastix_int_t M = pastix_imax(M2, M1);
    pastix_int_t N = pastix_imax(N2, N1);
    pastix_int_t minMN_1 = pastix_imin(M, rank);
    pastix_int_t minMN_2 = pastix_imin(N, rank);
    pastix_int_t i;
    pastix_complex64_t *u1u2, *v1v2, *R, *u, *v;
    pastix_complex64_t *tau1, *tau2;
    pastix_int_t ret;

    /* SVD entry parameters */
    double *s, *superb;
    pastix_complex64_t *tmp;
    pastix_int_t new_rank;
    Clock timer;

    clockStart(timer);
    assert(M2 == M);
    assert(N2 == N);

    /* Unused parameters right now */
    if (M1+offx > M2 || N1+offy > N2){
        errorPrint("Dimensions are not correct");
        return -1;
    }

    /* Rank is too high for u1u2 */
    if (minMN_1 == M){
        return -1;
    }

    /* Rank is too high for v1v2 */
    if (minMN_2 == N){
        return -1;
    }

    /**
     * Concatenate U2 and U1 in u1u2
     *  [ u2  0  ]
     *  [ u2  u1 ]
     *  [ u2  0  ]
     */
    u1u2 = malloc( M * rank * sizeof(pastix_complex64_t));
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, r2,
                         u2, ldu2, u1u2, M );

    tmp = u1u2 + r2 * M;
    if (u1 == NULL) {
        if (M1 != M2) {
            /* Set to 0 */
            memset(tmp, 0, M * r1 * sizeof(pastix_complex64_t));

            /* Set diagonal */
            tmp += offx;
            for (i=0; i<r1; i++, tmp += M+1) {
                *tmp = 1.;
            }
        }
        else {
            assert( offx == 0 );
            LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, r1,
                                 0., 1., tmp + offx, M );
        }
    }
    else{
        if (M1 != M) {
            memset(tmp, 0, M * r1 * sizeof(pastix_complex64_t));
        }
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M1, r1,
                             u1, ldu1, tmp + offx, M );
    }

    /**
     * Perform QR factorization on u1u2 = (Q1 R1)
     */
    tau1 = malloc( minMN_1 * sizeof(pastix_complex64_t));
    ret = LAPACKE_zgeqrf( LAPACK_COL_MAJOR, M, rank,
                          u1u2, M, tau1 );

    /**
     * Concatenate V2 and V1 in v1v2
     *  [ v2^h v2^h v2^h ]
     *  [ 0    v1^h 0    ]
     */
    v1v2 = malloc( N * rank * sizeof(pastix_complex64_t));
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', r2, N,
                         v2, ldv2, v1v2, rank );

    tmp = v1v2 + r2;
    if (v1 == NULL) {
        if (N1 != N2) {
            /* Set to 0 */
            LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', r1, N,
                                 0., 0., tmp, rank );

            /* Set diagonal */
            tmp += offy * rank;
            for (i=0; i<r1; i++, tmp += rank+1) {
                *tmp = alpha;
            }
        }
        else {
            assert( offy == 0 );
            LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', r1, N,
                                 0., alpha, tmp + offy * rank, rank );
        }
    }
    else{
        if (N1 != N) {
            LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', r1, N,
                                 0., 0., tmp, rank );
        }
        core_zgeadd( PastixNoTrans, r1, N1,
                     alpha, v1,                ldv1,
                        0., tmp + offy * rank, rank );
    }

    /**
     * Perform LQ factorization on v1v2 = (L2 Q2)
     */
    tau2 = malloc( minMN_2 * sizeof(pastix_complex64_t));
    ret = LAPACKE_zgelqf( LAPACK_COL_MAJOR, rank, N,
                          v1v2, rank, tau2 );

    /**
     * Compute R = alpha R1 L2
     */
    R = malloc( 3 * rank * rank * sizeof(pastix_complex64_t));
    u = R + rank * rank;
    v = u + rank * rank;

    memset(R, 0, rank * rank * sizeof(pastix_complex64_t));

    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'U', rank, rank,
                         u1u2, M, R, rank );

    cblas_ztrmm(CblasColMajor,
                CblasRight, CblasLower,
                CblasNoTrans, CblasNonUnit,
                rank, rank, CBLAS_SADDR(zone),
                v1v2, rank, R, rank);

    /**
     * Compute svd(R) = u \sigma v^t
     */
    s = malloc( rank * sizeof(double));
    superb = malloc( rank * sizeof(double));

    ret = LAPACKE_zgesvd( CblasColMajor, 'S', 'S',
                          rank, rank, R, rank,
                          s, u, rank, v, rank, superb );

    assert(ret == 0);
    if (ret != 0) {
        errorPrint("LAPACKE_zgesvd FAILED");
        EXIT(MOD_SOPALIN, INTERNAL_ERR);
    }

    new_rank = rank;
    for (i=0; i<rank; i++){
        if (s[i] < tol){ /// s[0] < tolerancep || s[i] < 0.5*tolerance){
            new_rank = i+1;
            break;
        }
    }

    /**
     * Let's now compute the final U = Q1 ([u] \sigma)
     *                                     [0]
     */
    /* First, let's scale the small u and compute the new rank */
    tmp = u;
    for (i=0; i<rank; i++, tmp+=rank){
        cblas_zscal(rank, CBLAS_SADDR(s[i]), tmp, 1);
    }
    if (new_rank == 0) {
        fprintf(stderr, "Rank null\n");
        /* errorPrint("Rank null\n", */
        /*            EXIT(MOD_SOPALIN, INTERNAL_ERR); */
    }

    /* Copy u at the top of u2 */
    tmp = u2;
    for (i=0; i<rank; i++, tmp+=ldu2, u+=rank) {
        memcpy(tmp, u,              rank * sizeof(pastix_complex64_t));
        memset(tmp + rank, 0, (M - rank) * sizeof(pastix_complex64_t));
    }

    ret = LAPACKE_zunmqr(LAPACK_COL_MAJOR, 'L', 'N',
                         M, rank, minMN_1,
                         u1u2, M, tau1,
                         u2,   ldu2);

    /**
     * And the final V^T = [v^t 0 ] Q2
     */
    assert( ldv2 >= new_rank );
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', rank, rank,
                         v, rank, v2, ldv2 );
    LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', rank, N-rank,
                         0., 0., v2 + ldv2 * rank, ldv2 );

    ret = LAPACKE_zunmlq(LAPACK_COL_MAJOR, 'R', 'N',
                         rank, N, minMN_2,
                         v1v2, rank, tau2,
                         v2, ldv2);

    free(u1u2);
    free(v1v2);
    free(s);
    free(superb);
    free(tau1);
    free(tau2);
    free(R);

    clockStop(timer);
    time_recomp += clockVal(timer);
    return new_rank;
}

pastix_int_t
core_z_add_LR3( double tol, pastix_complex64_t alpha,
                pastix_int_t M1, pastix_int_t N1, pastix_int_t r1,
                const pastix_complex64_t *u1, pastix_int_t ldu1,
                const pastix_complex64_t *v1, pastix_int_t ldv1,
                pastix_int_t M2, pastix_int_t N2, pastix_int_t r2,
                pastix_complex64_t *u2, pastix_int_t ldu2,
                pastix_complex64_t *v2, pastix_int_t ldv2,
                pastix_int_t offx, pastix_int_t offy)
{
    pastix_int_t rank  = r1 + r2;
    pastix_int_t M = pastix_imax(M1, M2);
    pastix_int_t N = pastix_imax(N1, N2);
    pastix_int_t minMN_u = pastix_imin(M, rank);
    pastix_int_t minMN_v = pastix_imin(N, rank);
    pastix_int_t i, ret, new_rank;
    pastix_complex64_t *u2u1, *v2v1, *R, *u, *v;
    pastix_complex64_t *tauU, *tauV, *tmp;

    /* SVD entry parameters */
    pastix_int_t j;
    pastix_complex64_t *Ru, *Rv;
    double *s, *superb;

    Clock timer;
    clockStart(timer);

    assert(M == M2);
    assert(N == N2);
    if (M1+offx > M2 || N1+offy > N2){
        errorPrint("Dimensions are not correct");
        return -1;
    }

    /* Rank is too high for u2u1 */
    if (minMN_u == M){
        return -1;
    }

    /* Rank is too high for v2v1 */
    if (minMN_v == N){
        return -1;
    }

    /* Now that we restrain the space of U, let's make sure we don't have a rank
     * larger than the space available
     * To be changed in a future version in add_LRv since this one will disappear */
    if (rank > pastix_imin(ldu2, ldv2) ) {
        return -1;
    }

    /**
     * Concatenate U2 and U1 in u1u2
     *  [ u2  0  ]
     *  [ u2  u1 ]
     *  [ u2  0  ]
     */
    u2u1 = malloc( M * rank * sizeof(pastix_complex64_t));
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, r2,
                         u2, ldu2, u2u1, M );

    tmp = u2u1 + r2 * M;
    if (u1 == NULL) {
        if (M1 != M2) {
            /* Set to 0 */
            memset(tmp, 0, M * r1 * sizeof(pastix_complex64_t));

            /* Set diagonal */
            tmp += offx;
            for (i=0; i<r1; i++, tmp += M+1) {
                *tmp = 1.;
            }
        }
        else {
            assert( offx == 0 );
            LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, r1,
                                 0., 1., tmp + offx, M );
        }
    }
    else{
        if (M1 != M) {
            memset(tmp, 0, M * r1 * sizeof(pastix_complex64_t));
        }
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M1, r1,
                             u1, ldu1, tmp + offx, M );
    }

    /**
     * Perform QR factorization on u2u1 = (Q1 R1)
     */
    tauU = malloc( minMN_u * sizeof(pastix_complex64_t));
    ret = LAPACKE_zgeqrf( LAPACK_COL_MAJOR, M, rank,
                          u2u1, M, tauU );


    /**
     * Concatenate V2 and V1 in v1v2
     *  [ v2  0  ]
     *  [ v2  v1 ]
     *  [ v2  0  ]
     */
    v2v1 = malloc( N * rank * sizeof(pastix_complex64_t));
    core_zgetro( r2, N, v2, ldv2, v2v1, N );

    tmp = v2v1 + r2 * N;
    if (v1 == NULL) {
        if (N1 != N2) {
            /* Set to 0 */
            memset(tmp, 0, N * r1 * sizeof(pastix_complex64_t));

            /* Set diagonal */
            tmp += offy;
            for (i=0; i<r1; i++, tmp += N+1) {
                *tmp = alpha;
            }
        }
        else {
            assert( offy == 0 );
            LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', N, r1,
                                 0., alpha, tmp + offy, N );
        }
    }
    else{
        if (N1 != N) {
            memset(tmp, 0, N * r1 * sizeof(pastix_complex64_t));
        }
        core_zgeadd( PastixConjTrans, N1, r1,
                     alpha, v1,         ldv1,
                     0.,    tmp + offy, N );
    }

    /**
     * Perform QR factorization on v2v1 = (Q1 R1)
     */
    tauV = malloc( minMN_v * sizeof(pastix_complex64_t));
    ret = LAPACKE_zgeqrf( LAPACK_COL_MAJOR, N, rank,
                          v2v1, N, tauV );

    Ru = malloc(rank * rank * sizeof(pastix_complex64_t));
    Rv = malloc(rank * rank * sizeof(pastix_complex64_t));
    R  = malloc(rank * rank * sizeof(pastix_complex64_t));
    memset(Ru, 0, rank * rank * sizeof(pastix_complex64_t));
    memset(Rv, 0, rank * rank * sizeof(pastix_complex64_t));

    for (i=0; i<rank; i++){
        memcpy(Ru + rank * i, u2u1 + M * i, (i+1) * sizeof(pastix_complex64_t));
        memcpy(Rv + rank * i, v2v1 + N * i, (i+1) * sizeof(pastix_complex64_t));
    }

    /* Compute Ru Rv^T */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                rank, rank, rank,
                CBLAS_SADDR(zone),  Ru, rank,
                                    Rv, rank,
                CBLAS_SADDR(zzero), R,  rank);

    s = malloc( rank * sizeof(double));
    u = malloc( rank * rank * sizeof(pastix_complex64_t));
    v = malloc( rank * rank * sizeof(pastix_complex64_t));
    superb = malloc( rank * sizeof(double));

    ret = LAPACKE_zgesvd( CblasColMajor, 'S', 'S',
                          rank, rank, R, rank,
                          s, u, rank, v, rank, superb );

    assert(ret == 0);
    if (ret != 0){
        errorPrint("LAPACKE_zgesvd FAILED");
    }

    new_rank = rank;
    for (i=0; i<rank; i++){
        if (s[i] < tol){ /// s[0] < tolerancep || s[i] < 0.5*tolerance){
            new_rank = i + 1;
            break;
        }
    }

    /* Scal u as before to take into account singular values */
    for (i=0; i<rank; i++){
        cblas_zscal(rank, CBLAS_SADDR(s[i]), u + rank * i, 1);
    }
    for (i=0; i<rank; i++){
        memcpy(u2 + M * i, u + rank * i, rank * sizeof(pastix_complex64_t));
        memset(u2 + M * i + rank, 0, (M2 - rank) * sizeof(pastix_complex64_t));
    }

    /* We need non-transposed version of v */
    pastix_complex64_t *v3;
    v3 = malloc( N * rank * sizeof(pastix_complex64_t) );

    for (i=0; i<rank; i++){
        for (j=0; j<rank; j++){
            v3[N * j + i] = v[rank * i + j];
        }
    }
    for (j=0; j<rank; j++){
        memset(v3 + N * j + rank, 0, (N - rank) * sizeof(pastix_complex64_t));
    }

    ret = LAPACKE_zunmqr(LAPACK_COL_MAJOR, 'L', 'N',
                         M, rank, minMN_u,
                         u2u1, M, tauU,
                         u2, M2);

    /* TODO: checker if 'R' 'T' can be realized */
    ret = LAPACKE_zunmqr(LAPACK_COL_MAJOR, 'L', 'N',
                         N, rank, minMN_v,
                         v2v1, N, tauV,
                         v3, N2);


    for (i=0; i<rank; i++){
        for (j=0; j<N2; j++){
            v2[N2 * j + i] = v3[N2 * i + j];
        }
    }

    free(u2u1);
    free(v2v1);
    free(u);
    free(s);
    free(v);
    free(superb);
    free(tauU);
    free(tauV);
    free(R);
    free(Ru);
    free(Rv);
    free(v3);

    clockStop(timer);
    time_recomp += clockVal(timer);
    return new_rank;
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
 *******************************************************************************
 */
pastix_int_t
core_z_add_LR(pastix_complex64_t *u1,
              pastix_complex64_t *v1,
              pastix_int_t dim_u1,
              pastix_int_t dim_v1,
              pastix_int_t rank_1,
              pastix_int_t ld_u1,
              pastix_int_t ld_v1,
              const pastix_complex64_t *u2,
              const pastix_complex64_t *v2,
              pastix_int_t trans,
              pastix_int_t dim_u2,
              pastix_int_t dim_v2,
              pastix_int_t rank_2,
              pastix_int_t ld_u2,
              pastix_int_t ld_v2,
              pastix_int_t x2,
              pastix_int_t y2)
{
    pastix_complex64_t alpha;

    /* Unused parameters right now */
    (void)ld_u1;
    (void)ld_v1;

    Clock timer;
    clockStart(timer);

    pastix_int_t dim_u = pastix_imax(dim_u1, dim_u2);
    pastix_int_t dim_v = pastix_imax(dim_v1, dim_v2);
    pastix_int_t rank  = rank_1 + rank_2;
    pastix_int_t i, j;

    if (dim_u2 > dim_u1 || dim_v2 > dim_v1){
        errorPrint("Dimensions are not correct");
    }

    pastix_int_t minMN_1 = pastix_imin(dim_u, rank);
    /* Rank is too high for u1u2 */
    if (minMN_1 == dim_u){
        return -1;
    }

    pastix_int_t minMN_2 = pastix_imin(dim_v, rank);
    /* Rank is too high for v1v2 */
    if (minMN_2 == dim_v){
        return -1;
    }

    /* Now that we restrain the space of U, let's make sure we don't have a rank
     * larger than the space available
     * To be changed in a future version in add_LR2 since this one will disappear */
    if (rank > pastix_imin(ld_u1, ld_v1) ) {
        return -1;
    }

    pastix_complex64_t *u1u2, *v1v2;
    pastix_int_t ret;
    pastix_complex64_t *tau1;
    pastix_complex64_t *tau2;

    /* Matrices for computing SVD of R1 R2^T */
    pastix_complex64_t *R1, *R2, *R;

    /* SVD entry parameters */
    double *s, *superb;
    pastix_complex64_t *u, *v;

    pastix_int_t new_rank;

    u1u2 = malloc( dim_u * rank * sizeof(pastix_complex64_t));
    v1v2 = malloc( dim_v * rank * sizeof(pastix_complex64_t));

    /* Complete with zeroes if adding a small block in a larger structure */
    if (dim_u1 != dim_u2)
        memset(u1u2 + dim_u1 * rank_1, 0, dim_u1 * rank_2 * sizeof(pastix_complex64_t));
    if (dim_v1 != dim_v2)
        memset(v1v2, 0, dim_v * rank * sizeof(pastix_complex64_t));

    memcpy(u1u2, u1, dim_u1 * rank_1 * sizeof(pastix_complex64_t));

    /* For identity matrix */
    if (u2 == NULL){
        for (i=0; i<rank_2; i++){
            u1u2[dim_u * (rank_1 + i) + x2 + i] = 1;
        }
    }
    else{
        for (i=0; i<rank_2; i++){
            memcpy(u1u2 + dim_u * (rank_1 + i) + x2, u2 + i * ld_u2, dim_u2 * sizeof(pastix_complex64_t));
        }
    }

    for (i=0; i<dim_v1; i++){
        for (j=0; j<rank_1; j++){
            v1v2[dim_v * j + i] = v1[i * dim_v1 + j];
        }
    }

    /* For identity matrix */
    if (v2 == NULL){
        for (j=0; j<rank_2; j++){
            v1v2[dim_v * (rank_1 + j) + j + y2] = -1;
        }
    }
    /* WARNING: minus because of the extend add */
    else{
        for (i=0; i<dim_v2; i++){
            for (j=0; j<rank_2; j++){
                if (trans == CblasNoTrans)
                    v1v2[dim_v * (rank_1 + j) + i + y2] = -v2[i * ld_v2 + j];
                else
                    v1v2[dim_v * (rank_1 + j) + i + y2] = -v2[j * ld_v2 + i];
            }
        }
    }

    tau1 = malloc( minMN_1 * sizeof(pastix_complex64_t));
    ret = LAPACKE_zgeqrf( CblasColMajor, dim_u, rank,
                          u1u2, dim_u, tau1 );


    tau2 = malloc( minMN_2 * sizeof(pastix_complex64_t));
    ret = LAPACKE_zgeqrf( CblasColMajor, dim_v, rank,
                          v1v2, dim_v, tau2 );

    R1 = malloc(rank * rank * sizeof(pastix_complex64_t));
    R2 = malloc(rank * rank * sizeof(pastix_complex64_t));
    R  = malloc(rank * rank * sizeof(pastix_complex64_t));
    memset(R1, 0, rank * rank * sizeof(pastix_complex64_t));
    memset(R2, 0, rank * rank * sizeof(pastix_complex64_t));

    for (i=0; i<rank; i++){
        memcpy(R1 + rank * i, u1u2 + dim_u * i, (i+1) * sizeof(pastix_complex64_t));
        memcpy(R2 + rank * i, v1v2 + dim_v * i, (i+1) * sizeof(pastix_complex64_t));
    }

    /* Compute R1 R2^T */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                rank, rank, rank,
                CBLAS_SADDR(zone),  R1, rank,
                                    R2, rank,
                CBLAS_SADDR(zzero), R,  rank);

    s = malloc( rank * sizeof(double));
    u = malloc( rank * rank * sizeof(pastix_complex64_t));
    v = malloc( rank * rank * sizeof(pastix_complex64_t));
    superb = malloc( rank * sizeof(double));

    ret = LAPACKE_zgesvd( CblasColMajor, 'S', 'S',
                          rank, rank, R, rank,
                          s, u, rank, v, rank, superb );

    assert(ret == 0);
    if (ret != 0){
        errorPrint("LAPACKE_zgesvd FAILED");
    }

    new_rank = rank;
    char *tol             = getenv("TOLERANCE");
    double tolerance      = atof(tol);

    for (i=0; i<rank-1; i++){
        if (s[i] < tolerance){ /// s[0] < tolerancep || s[i] < 0.5*tolerance){
            new_rank = i + 1;
            break;
        }
    }

    /* Scal u as before to take into account singular values */
    for (i=0; i<rank; i++){
        alpha = s[i];
        cblas_zscal(rank, CBLAS_SADDR(alpha), u + rank * i, 1);
    }
    for (i=0; i<rank; i++){
        memcpy(u1 + dim_u * i, u + rank * i, rank * sizeof(pastix_complex64_t));
        memset(u1 + dim_u * i + rank, 0, (dim_u1 - rank) * sizeof(pastix_complex64_t));
    }

    /* We need non-transposed version of v */
    pastix_complex64_t *v3;
    v3 = malloc( dim_v * rank * sizeof(pastix_complex64_t) );

    for (i=0; i<rank; i++){
        for (j=0; j<rank; j++){
            v3[dim_v * j + i] = v[rank * i + j];
        }
    }
    for (j=0; j<rank; j++){
        memset(v3 + dim_v * j + rank, 0, (dim_v - rank) * sizeof(pastix_complex64_t));
    }

    ret = LAPACKE_zunmqr(LAPACK_COL_MAJOR, 'L', 'N',
                         dim_u, rank, minMN_1,
                         u1u2, dim_u, tau1,
                         u1, dim_u1);

    /* TODO: checker if 'R' 'T' can be realized */
    ret = LAPACKE_zunmqr(LAPACK_COL_MAJOR, 'L', 'N',
                         dim_v, rank, minMN_2,
                         v1v2, dim_v, tau2,
                         v3, dim_v1);


    for (i=0; i<rank; i++){
        for (j=0; j<dim_v1; j++){
            v1[dim_v1 * j + i] = v3[dim_v1 * i + j];
        }
    }

    free(u1u2);
    free(v1v2);
    free(u);
    free(s);
    free(v);
    free(superb);
    free(tau1);
    free(tau2);
    free(R);
    free(R1);
    free(R2);
    free(v3);

    clockStop(timer);
    time_recomp += clockVal(timer);
    return new_rank;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_z_lr2dense - If the u v^T matrix is still used, uncompress it into A
 *
 *
 *******************************************************************************/
void core_z_lr2dense(SolverBlok *blok,
                     pastix_complex64_t *A,
                     pastix_int_t stride,
                     pastix_int_t width,
                     pastix_int_t side){

    switch (side){
    case L_side:
        if (blok->rankL != -1){
            pastix_int_t dimb     = blok->lrownum - blok->frownum + 1;
            pastix_complex64_t *u = blok->coefL_u_LR;
            pastix_complex64_t *v = blok->coefL_v_LR;

            core_z_uncompress_LR(dimb, width, blok->rankL,
                                 u, dimb,
                                 v, width,
                                 A, stride);
        }
        break;
    case U_side:
        if (blok->rankU != -1){
            pastix_int_t dimb     = blok->lrownum - blok->frownum + 1;
            pastix_complex64_t *u = blok->coefU_u_LR;
            pastix_complex64_t *v = blok->coefU_v_LR;

            core_z_uncompress_LR(dimb, width, blok->rankU,
                                 u, dimb,
                                 v, width,
                                 A, stride);
        }
        break;
    default:
        errorPrint("Wrong operation");
        break;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zproduct_lr2dense - Compute work = A1 * A2
 * A1, A2 maybe dense or LR
 * work has to be dense
 *
 *
 *******************************************************************************/
void core_zproduct_lr2dense(SolverBlok *blok1,
                            pastix_complex64_t *A1,
                            pastix_int_t stride1,
                            pastix_int_t width1,
                            pastix_int_t side1,
                            SolverBlok *blok2,
                            pastix_complex64_t *A2,
                            pastix_int_t stride2,
                            pastix_int_t width2,
                            pastix_int_t side2,
                            pastix_complex64_t *work,
                            pastix_int_t ldwork){

    /* TODO: enhance if A1 and/or A2 are really LR */
    assert(width1 == width2);

    pastix_int_t dimb = blok1->lrownum - blok1->frownum + 1;
    pastix_int_t dimj = blok2->lrownum - blok2->frownum + 1;
    pastix_int_t dima = width1;

    core_z_lr2dense(blok1, A1, stride1,
                    width1, side1);

    core_z_lr2dense(blok2, A2, stride2,
                    width2, side2);

    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                 dimb, dimj, dima,
                 CBLAS_SADDR(zone),  A1, stride1,
                                     A2, stride2,
                 CBLAS_SADDR(zzero), work, ldwork  );
}

pastix_int_t
core_zlradd(double tol, pastix_complex64_t alpha,
            pastix_int_t M1, pastix_int_t N1, pastix_int_t r1,
            const pastix_complex64_t *u1, pastix_int_t ldu1,
            const pastix_complex64_t *v1, pastix_int_t ldv1,
            pastix_int_t M2, pastix_int_t N2, pastix_int_t r2,
            pastix_complex64_t *u2, pastix_int_t ldu2,
            pastix_complex64_t *v2, pastix_int_t ldv2,
            pastix_int_t offx, pastix_int_t offy)
{
    //#define OLD_EXTENDADD
#if defined(OLD_EXTENDADD)
    (void)alpha; (void)tol;
    return core_z_add_LR( u2, v2,
                          M2, N2, r2, ldu2, ldv2,
                          u1, v1, CblasNoTrans,
                          M1, N1, r1, ldu1, ldv1,
                          offx, offy);
#else
    return core_z_add_LR2( tol, alpha,
                           M1, N1, r1, u1, ldu1, v1, ldv1,
                           M2, N2, r2, u2, ldu2, v2, ldv2,
                           offx, offy);
#endif
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zproduct_lr2lr - Update uF vF^T with A1 * A2
 * A1, A2 maybe dense or LR
 * uf vF^T is necessary LR
 *
 *
 *******************************************************************************/
void core_zproduct_lr2lr(SolverBlok *blok1,
                         pastix_complex64_t *A1,
                         pastix_int_t stride1,
                         pastix_int_t width1,
                         pastix_int_t side1,
                         SolverBlok *blok2,
                         pastix_complex64_t *A2,
                         pastix_int_t stride2,
                         pastix_int_t width2,
                         pastix_int_t side2,
                         SolverBlok *blok3,
                         pastix_complex64_t *A3,
                         pastix_complex64_t *A3_contrib,
                         pastix_int_t stride3,
                         pastix_int_t width3,
                         pastix_int_t side3,
                         pastix_int_t offset,
                         pastix_complex64_t *work){

    pastix_int_t dimb = blok1->lrownum - blok1->frownum + 1;
    pastix_int_t dimj = blok2->lrownum - blok2->frownum + 1;
    pastix_int_t dima = width1;

    pastix_complex64_t *uL, *vL, *uU, *vU, *uF, *vF;
    //pastix_complex64_t *Au, *Av, *Bu, *Bv, *Cu, *Cv;
    pastix_int_t rankL, rankU, rankF;
    //pastix_int_t rkA, rkB, rkC;

    char *tolerance = getenv("TOLERANCE");
    double tol = atof(tolerance);

    rankL = -1;
    rankU = -1;
    rankF = -1;

    uL = NULL;
    vL = NULL;
    uU = NULL;
    vU = NULL;
    uF = NULL;
    vF = NULL;

    /* Operation L * U contributes to L */
    if (side1 == L_side && side2 == U_side && side3 == L_side){
        rankL = blok1->rankL;
        rankU = blok2->rankU;
        rankF = blok3->rankL;

        uL = blok1->coefL_u_LR;
        vL = blok1->coefL_v_LR;
        uU = blok2->coefU_u_LR;
        vU = blok2->coefU_v_LR;
        uF = blok3->coefL_u_LR;
        vF = blok3->coefL_v_LR;
    }
    else if (side1 == U_side && side2 == L_side && side3 == U_side){
        rankL = blok1->rankU;
        rankU = blok2->rankL;
        rankF = blok3->rankU;

        uL = blok1->coefU_u_LR;
        vL = blok1->coefU_v_LR;
        uU = blok2->coefL_u_LR;
        vU = blok2->coefL_v_LR;
        uF = blok3->coefU_u_LR;
        vF = blok3->coefU_v_LR;
    }
    else{
        errorPrint("Operation not supported yet");
    }

    /* Au = uL; */
    /* Av = vL; */
    /* Bu = uU; */
    /* Bv = vU; */
    /* Cu = uF; */
    /* Cv = vF; */
    /* rkA = rankL; */
    /* rkB = rankU; */
    /* rkC = rankF; */

    pastix_int_t ret   = -1;
    pastix_int_t sizeF = blok3->lrownum - blok3->frownum + 1;

    /* Perform LR * LR */
    if (rankU != -1 && rankL != -1){

        pastix_complex64_t *tmp;
        pastix_int_t stride_tmp = pastix_imax(dimj, dimb);
        tmp = malloc(dimj * stride_tmp * sizeof(pastix_complex64_t *));

        /**
         * Let's compute A * B = Au Av^h (Bu Bv^h)' in two steps
         *    1) w = Av^h Bv^h' that is normally the smaller of the three products
         *    2) If M < N, Au * w, or w * Bu' otherwise
         */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                     rankL, rankU, dima,
                     CBLAS_SADDR(zone),  vL,   dima,
                                         vU,   dima,
                     CBLAS_SADDR(zzero), work, dimb );

        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                     rankL, dimj, rankU,
                     CBLAS_SADDR(zone),  work, dimb,
                                         uU,   dimj,
                     CBLAS_SADDR(zzero), tmp,  stride_tmp );

        ret = core_zlradd( tol, -1,
                           /* A*B */
                           dimb, dimj, rankL,
                           uL, dimb, tmp, stride_tmp,
                           /* C */
                           sizeF, width3, rankF,
                           uF, sizeF, vF, width3,
                           /* offset */
                           blok1->frownum - blok3->frownum,
                           blok2->frownum - offset);

        /* TODO: enhance, but requires memory... */
        /* Form A1 A2 for dense contribution */
        if (ret == -1){
            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         dimb, dimj, rankL,
                         CBLAS_SADDR(zone),  uL,   dimb,
                                             tmp,  stride_tmp,
                         CBLAS_SADDR(zzero), work, dimb  );
        }

        free(tmp);
    }

    /* Perform dense matrix */
    else if (rankU == -1 && rankL == -1){

        if (dima < dimb && dima < dimj){
            /* Add matrix of rank dima */
            ret = core_z_add_LR(uF, vF,
                                sizeF, width3, rankF, sizeF, width3,
                                A1, A2, CblasTrans,
                                dimb, dimj, dima,
                                stride1, stride2,
                                blok1->frownum - blok3->frownum,
                                blok2->frownum - offset);

            if (ret == -1){
                /* A1 A2 matrix-matrix product was not form */
                cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                             dimb, dimj, dima,
                             CBLAS_SADDR(zone),  A1,  stride1,
                             A2,  stride2,
                             CBLAS_SADDR(zzero), work, dimb );
            }
        }
        else{
            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                         dimb, dimj, dima,
                         CBLAS_SADDR(zone),  A1,  stride1,
                         A2,  stride2,
                         CBLAS_SADDR(zzero), work, dimb );

            if (dimb < dimj){
                /* Add matrix of rank dimb */
                ret = core_zlradd( tol, -1,
                                   /* A*B */
                                   dimb, dimj, dimb,
                                   NULL, dimb, work, dimb,
                                   /* C */
                                   sizeF, width3, rankF,
                                   uF, sizeF, vF, width3,
                                   /* offset */
                                   blok1->frownum - blok3->frownum,
                                   blok2->frownum - offset);
            }
            else{
                /* Add matrix of rank dimj */
                ret = core_zlradd( tol, -1,
                                   /* A*B */
                                   dimb, dimj, dimj,
                                   work, dimb, NULL, dimj,
                                   /* C */
                                   sizeF, width3, rankF,
                                   uF, sizeF, vF, width3,
                                   /* offset */
                                   blok1->frownum - blok3->frownum,
                                   blok2->frownum - offset);
            }
        }
    }

    /* Perform dense * LR */
    else if (rankU == -1 && rankL != -1){
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                     rankL, dimj, dima,
                     CBLAS_SADDR(zone),  vL, dima,
                     A2, stride2,
                     CBLAS_SADDR(zzero), work, dimb  );

        ret = core_zlradd( tol, -1,
                           /* A*B */
                           dimb, dimj, rankL,
                           uL, dimb, work, dimb,
                           /* C */
                           sizeF, width3, rankF,
                           uF, sizeF, vF, width3,
                           /* offset */
                           blok1->frownum - blok3->frownum,
                           blok2->frownum - offset);

        /* TODO: enhance, but requires memory... */
        /* Form A1 A2 for dense contribution */
        if (ret == -1){
            core_z_lr2dense(blok1, A1, stride1,
                            width1, side1);

            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                         dimb, dimj, dima,
                         CBLAS_SADDR(zone),  A1,  stride1,
                         A2,  stride2,
                         CBLAS_SADDR(zzero), work, dimb );
        }
    }

    /* Perform LR * dense */
    else if (rankU != -1 && rankL == -1){
        /* TODO: Change to Al * (Av^t * B) or you compute dense * LR ????? */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                     dimb, rankU, dima,
                     CBLAS_SADDR(zone),  A1, stride1,
                                         vU, dima,
                     CBLAS_SADDR(zzero), work, dimb  );

        ret = core_z_add_LR(uF, vF,
                            sizeF, width3, rankF, sizeF, width3,
                            work, uU, CblasTrans,
                            dimb, dimj, rankU, dimb, dimj,
                            blok1->frownum - blok3->frownum,
                            blok2->frownum - offset);

        /* TODO: enhance, but requires memory... */
        /* Form A1 A2 for dense contribution */
        if (ret == -1){
            core_z_lr2dense(blok2, A2, stride2,
                            width2, side2);

            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                         dimb, dimj, dima,
                         CBLAS_SADDR(zone),  A1,  stride1,
                         A2,  stride2,
                         CBLAS_SADDR(zzero), work, dimb );
        }
    }


    /* Uncompress the block and perform dense update if rank became too large */
    if (ret == -1){
        core_z_lr2dense(blok3, A3, stride3,
                        width3, side3);

        core_zgeadd( CblasNoTrans, dimb, dimj,
                     -1.0, work, dimb,
                      1.0, A3_contrib, stride3 );
    }
    if (side3 == L_side){
        blok3->rankL = ret;
    }
    else if (side3 == U_side){
        blok3->rankU = ret;
    }
}


/* We add current zone by modifying the beta parameter in gemm */
void core_z_lr2dense_special(SolverBlok *blok,
                             pastix_complex64_t *A,
                             pastix_int_t stride,
                             pastix_int_t width,
                             pastix_int_t side){

    switch (side){
    case L_side:
        if (blok->rankL != -1){
            pastix_int_t dimb     = blok->lrownum - blok->frownum + 1;
            pastix_complex64_t *u = blok->coefL_u_LR;
            pastix_complex64_t *v = blok->coefL_v_LR;
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        dimb, width, blok->rankL,
                        CBLAS_SADDR(zone),  u,  dimb,
                                            v,  width,
                        CBLAS_SADDR(mzone), A, stride);
        }
        break;
    case U_side:
        if (blok->rankU != -1){
            pastix_int_t dimb     = blok->lrownum - blok->frownum + 1;
            pastix_complex64_t *u = blok->coefU_u_LR;
            pastix_complex64_t *v = blok->coefU_v_LR;
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        dimb, width, blok->rankU,
                        CBLAS_SADDR(zone),  u,  dimb,
                                            v,  width,
                        CBLAS_SADDR(mzone), A, stride);
        }
        break;
    default:
        errorPrint("Wrong operation");
        break;
    }
}
