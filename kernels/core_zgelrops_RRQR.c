/**
 *
 * @file core_zgelrops_RRQR.c
 *
 * PaStiX low-rank kernel routines using Rank-revealing QR based on Lapack GEQP3.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Alfredo Buttari
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
static pastix_complex64_t mzone = -1.;
static pastix_complex64_t zone  =  1.;
static pastix_complex64_t zzero =  0.;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Compute a rank-reavealing QR factorization.
 *
 * This routine is originated from the LAPACK kernels zgeqp3/zlaqps and was
 * modified by A. Buttari for MUMPS-BLR.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          Number of rows of the matrix A.
 *
 * @param[in] n
 *          Number of columns of the matrix A.
 *
 * @param[in] A
 *          The matrix of dimension lda-by-n that need to be compressed
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1, m)
 *
 * @param[out] jpvt
 *          The array that describes the permutation of A
 *
 * @param[out] tau
 *          Contains scalar factors of the elementary reflectors for the matrix
 *          Q
 *
 * @param[in] work
 *          Workspace array
 *
 * @param[in] ldwork
 *          The leading dimension of work
 *
 * @param[in] rwork
 *          Workspace array used to store partial and exact column norms
 *
 * @param[in] tol
 *          The relative tolerance criteria. Computations are stopped when the
 *          norm of the residual matrix is lower than tol.
 *
 * @param[in] nb
 *          Blocking size for GEMM
 *
 * @param[in] maxrank
 *         Maximum rank of A. Computations are stopped when the rank exceeds
 *         maxrank
 *
 *******************************************************************************
 *
 * @return  This routine will return the rank of A (>=0) or an error (<0)
 *
 *******************************************************************************/
int
core_zrrqr( pastix_int_t m, pastix_int_t n,
            pastix_complex64_t *A, pastix_int_t lda,
            pastix_int_t *jpvt, pastix_complex64_t *tau,
            pastix_complex64_t *work, pastix_int_t ldwork,
            double *rwork,
            double tol, pastix_int_t nb, pastix_int_t maxrank)
{
    pastix_int_t minMN = pastix_imin(m, n);
    pastix_int_t ldf   = ldwork;

    pastix_int_t j, k, jb, itemp, lsticc, pvt;
    double temp, temp2;
    double machine_prec = sqrt(LAPACKE_dlamch('e'));
    pastix_complex64_t akk;

    pastix_complex64_t *auxv, *f;

    /* Partial (VN1) and exact (VN2) column norms */
    double *VN1, *VN2;

    /* Number or rows of A that have been factorized */
    pastix_int_t offset = 0;

    /* Rank */
    pastix_int_t rk = 0;

    if (m < 0)
        return -1;
    if (n < 0)
        return -2;
    if (lda < pastix_imax(1, m))
        return -4;
    if( ldwork < n)
        return -8;

    VN1 = rwork;
    VN2 = rwork + n;

    auxv = work;
    f    = work + ldwork;

    /* Initialize partial column norms. The first N elements of work */
    /* store the exact column norms. */
    /* TODO: call PLASMA/internal kernel */
    for (j=0; j<n; j++){
        VN1[j]  = cblas_dznrm2(m, A + j*lda, 1);
        VN2[j]  = VN1[j];
        jpvt[j] = j;
    }

    while(rk <= maxrank) {
        /* jb equivalent to kb in LAPACK xLAQPS: number of columns actually factorized */
        jb     = pastix_imin(nb, minMN-offset+1);
        lsticc = 0;

        /* column being factorized among jb */
        k = 0;

        while(k < jb && lsticc == 0){

            rk = offset+k;

            /* Rank is too large for compression */
            if (rk > maxrank){
                return rk;
            }

            pvt = rk + cblas_izamax(n-rk, VN1 + rk, 1);

            if (VN1[pvt] <= tol){
                double residual = 0.;
                pastix_int_t i;
                for (i=rk; i<n; i++){
                    residual += VN1[i]*VN1[i];
                }
                if (sqrt(residual) <= tol)
                    return rk;
            }

            /* Pivot is not within the current column: we swap */
            if (pvt != rk){
                assert( (pvt < n) && (rk < n) );
                cblas_zswap(m, A + pvt * lda, 1, A + rk * lda, 1);
                cblas_zswap(k, f + (pvt-offset), ldf, f + k, ldf);

                itemp     = jpvt[pvt];
                jpvt[pvt] = jpvt[rk];
                jpvt[rk]  = itemp;
                VN1[pvt]  = VN1[rk];
                VN2[pvt]  = VN2[rk];
            }

            /* Apply previous Householder reflectors to column K */
            /* A(RK:M,RK) := A(RK:M,RK) - A(RK:M,OFFSET+1:RK-1)*F(K,1:K-1)**H */
            if (k > 0){
#if defined(PRECISION_c) || defined(PRECISION_z)
                for (j=0; j<k; j++){
                    f[j * ldf + k] = conj(f[j * ldf + k]);
                }
#endif

                assert( (offset+1) < n );
                assert( (rk < n) && (rk < m) );
                cblas_zgemv(CblasColMajor, CblasNoTrans, m-rk, k, CBLAS_SADDR(mzone),
                            A + (offset) * lda + rk, lda,
                            f + k, ldf,
                            CBLAS_SADDR(zone), A + rk * lda + rk, 1);

#if defined(PRECISION_c) || defined(PRECISION_z)
                for (j=0; j<k; j++){
                    f[j * ldf + k] = conj(f[j * ldf + k]);
                }
#endif
            }

            /* Generate elementary reflector H(k). */
            if (rk < (m-1)){
                LAPACKE_zlarfg(m-rk, A + rk * lda + rk, A + rk * lda + (rk+1), 1, tau + rk);
            }
            else{
                LAPACKE_zlarfg(1, A + rk * lda + rk, A + rk * lda + rk, 1, tau + rk);
            }

            akk = A[rk * lda + rk];
            A[rk * lda + rk] = zone;

            /* Compute Kth column of F: */
            /* F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K). */
            if (rk < (n-1)){
                pastix_complex64_t alpha = tau[rk];
                cblas_zgemv(CblasColMajor, CblasConjTrans, m-rk, n-rk-1, CBLAS_SADDR(alpha),
                            A + (rk+1) * lda + rk, lda,
                            A + rk * lda + rk, 1,
                            CBLAS_SADDR(zzero), f + k * ldf + k + 1, 1);
            }

            /* Padding F(1:K,K) with zeros. */
            for (j=0; j<k; j++){
                f[k * ldf + j] = zzero;
            }

            /* Incremental updating of F: */
            /* F(1:N,K) := F(1:N-OFFSET,K) - tau(RK)*F(1:N,1:K-1)*A(RK:M,OFFSET+1:RK-1)**H*A(RK:M,RK). */
            if (k > 0){
                pastix_complex64_t alpha = -tau[rk];
                cblas_zgemv(CblasColMajor, CblasConjTrans, m-rk, k, CBLAS_SADDR(alpha),
                            A + (offset) * lda + rk, lda,
                            A + rk * lda + rk, 1,
                            CBLAS_SADDR(zzero), auxv, 1);

                cblas_zgemv(CblasColMajor, CblasNoTrans, n-offset, k, CBLAS_SADDR(zone),
                            f, ldf,
                            auxv, 1,
                            CBLAS_SADDR(zone), f + k * ldf, 1);
            }

            /* Update the current row of A: */
            /* A(RK,RK+1:N) := A(RK,RK+1:N) - A(RK,OFFSET+1:RK)*F(K+1:N,1:K)**H. */
            if (rk < (n-1)){

#if defined(PRECISION_c) || defined(PRECISION_z)
                cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                            1, n-rk-1, k+1,
                            CBLAS_SADDR(mzone), A + (offset) * lda + rk, lda,
                            f + (k+1), ldf,
                            CBLAS_SADDR(zone), A + (rk+1) * lda + rk, lda);
#else
                cblas_zgemv(CblasColMajor, CblasNoTrans, n-rk-1, k+1, CBLAS_SADDR(mzone),
                            f + (k+1), ldf,
                            A + (offset) * lda + rk, lda,
                            CBLAS_SADDR(zone), A + (rk+1) * lda + rk, lda);
#endif
            }

            /* Update partial column norms. */
            if (rk < (minMN-1)){
                for (j=rk+1; j<n; j++){
                    if (VN1[j] != 0.0){
                        /* NOTE: The following 4 lines follow from the analysis in */
                        /* Lapack Working Note 176. */
                        temp  = cabs( A[j * lda + rk] ) / VN1[j];
                        double temp3 = (1.0 + temp) * (1.0 - temp);
                        if (temp3 > 0.0){
                            temp = temp3;
                        }
                        else{
                            temp = 0.;
                        }
                        temp2 = temp * ( VN1[j] / VN2[j]) * ( VN1[j] / VN2[j]);
                        if (temp2 < machine_prec){
                            /* printf("LSTICC %ld\n", j); */
                            VN2[j] = lsticc;
                            lsticc = j;
                        }
                        else{
                            VN1[j] = VN1[j] * sqrt(temp);
                        }

                    }
                }
            }

            A[rk * lda + rk] = akk;

            k++;
        }

        /* Apply the block reflector to the rest of the matrix: */
        /* A(RK+1:M,RK+1:N) := A(RK+1:M,RK+1:N) - */
        /* A(RK+1:M,OFFSET+1:RK)*F(K+1:N-OFFSET,1:K)**H. */
        if (rk < (minMN-1)){
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                        m-rk-1, n-rk-1, k,
                        CBLAS_SADDR(mzone), A + (offset) * lda + rk + 1, lda,
                        f + k, ldf,
                        CBLAS_SADDR(zone), A + (rk+1) * lda + rk + 1, lda);
        }

        /* Recomputation of difficult columns. */
        while (lsticc > 0){
            itemp = (pastix_int_t) (VN2[lsticc]);
            assert(lsticc < n);
            VN1[lsticc] = cblas_dznrm2(m-rk-1, A + (lsticc) * lda + rk + 1, 1);

            /* NOTE: The computation of VN1( LSTICC ) relies on the fact that  */
            /* SNRM2 does not fail on vectors with norm below the value of */
            /* SQRT(DLAMCH('S'))  */
            VN2[lsticc] = VN1[lsticc];
            lsticc = itemp;
        }

        lsticc = 0;
        offset = rk+1;
    }

    return rk+1;
}

/**
 *******************************************************************************
 *
 * @brief Convert a full rank matrix in a low rank matrix, using RRQR.
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
void
core_zge2lr_RRQR( double tol, pastix_int_t m, pastix_int_t n,
                  const pastix_complex64_t *A, pastix_int_t lda,
                  pastix_lrblock_t *Alr )
{
    int                 ret;
    pastix_int_t        nb     = 32;
    pastix_int_t        ldwork = pastix_imax(m, n);
    pastix_complex64_t *work, *Acpy, *tau;
    double             *rwork;
    pastix_int_t       *jpvt;
    pastix_int_t        zsize, rsize;
    pastix_complex64_t *zwork;

    double norm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                       A, lda, NULL );

    /* work */
    zsize = (2 * nb + 1) * ldwork;
    /* Acpy */
    zsize += m * n;
    /* tau */
    zsize += n;

    /* rwork */
    rsize = 2* n;

    zwork = malloc( zsize * sizeof(pastix_complex64_t) + rsize * sizeof(double) );
    rwork = (double*)(zwork + zsize);

    work = zwork;
    Acpy = zwork + (2 * nb + 1) * ldwork;
    tau  = Acpy + m * n;

    MALLOC_INTERN( jpvt, n, pastix_int_t );

    /**
     * Allocate a temorary Low rank matrix
     */
    core_zlralloc( m, n, pastix_imin( m, n ), Alr );

    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                               A, lda, Acpy, m );
    assert(ret == 0);

    ret = core_zrrqr( m, n,
                      Acpy, m,
                      jpvt, tau,
                      work, ldwork,
                      rwork,
                      tol * norm, nb, pastix_imin(m,n) - 1 );

    /**
     * Resize the space used by the low rank matrix
     */
    core_zlrsze( 0, m, n, Alr, ret, -1 );

    /**
     * It was not interesting to compress, so we store the dense version in Alr
     */
    if ( ( Alr->rk == -1 ) ) {
        memFree_null(Alr->u);
        memFree_null(Alr->v);
        Alr->u     = malloc( m * n * sizeof(pastix_complex64_t) );
        Alr->rk    = -1;
        Alr->rkmax = m;

        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                                   A, lda, Alr->u, Alr->rkmax );
        assert(ret == 0);
    }
    /**
     * We compute Q/R obtained thanks to core_zrrqr
     */
    else if (Alr->rk != 0) {
        pastix_int_t i;

        /* Temporary space to permute Alr->v */
        pastix_complex64_t *permQR;
        pastix_complex64_t *U = Alr->u;
        pastix_complex64_t *V = Alr->v;

        MALLOC_INTERN( permQR, n * Alr->rk, pastix_complex64_t );
        memset( permQR, 0, n * Alr->rk * sizeof(pastix_complex64_t) );

        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'U', Alr->rk, n,
                                   Acpy, m,
                                   permQR, Alr->rk );
        assert(ret == 0);

        /* Permute V */
        for (i=0; i<n; i++){
            memcpy(V + jpvt[i] * Alr->rk,
                   permQR + i * Alr->rk,
                   Alr->rk * sizeof(pastix_complex64_t));
        }

        memFree_null( permQR );

        /* Compute Q factor on u */
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, Alr->rk,
                                   Acpy, m, U, m );
        assert(ret == 0);

        ret = LAPACKE_zungqr( LAPACK_COL_MAJOR, m, Alr->rk, Alr->rk,
                              U , m, tau );
        assert(ret == 0);
    }

    memFree_null( zwork );
    memFree_null( jpvt );
}

/**
 *******************************************************************************
 *
 * @brief Add two LR structures A=(-u1) v1^T and B=u2 v2^T into u2 v2^T
 *
 *    u2v2^T - u1v1^T = (u2 u1) (v2 v1)^T
 *    Orthogonalize (u2 u1) = (u2, u1 - u2(u2^T u1)) * (I u2^T u1)
 *                                                     (0    I   )
 *    Compute RRQR decomposition of (I u2^T u1) * (v2 v1)^T
 *                                  (0    I   )
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          The absolute tolerance criteria
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
 * @return  The new rank of u2 v2^T or -1 if ranks are too large for
 *          recompression
 *
 *******************************************************************************/
int
core_zrradd_RRQR( double tol, pastix_trans_t transA1, pastix_complex64_t alpha,
                  pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                  pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                  pastix_int_t offx, pastix_int_t offy)
{
    pastix_int_t rank, M, N, minV;
    pastix_int_t i, ret, new_rank;
    pastix_int_t ldau, ldav, ldbu, ldbv;
    pastix_complex64_t *u1u2, *v1v2, *u, *v;
    pastix_complex64_t *tmp, *zbuf, *tauV;
    size_t wzsize;
    double norm;

    /* RRQR parameters / workspace */
    pastix_int_t        nb = 32;
    pastix_int_t        ldwork;
    pastix_int_t       *jpvt;
    pastix_complex64_t *zwork;
    double             *rwork;

    rank = (A->rk == -1) ? pastix_imin(M1, N1) : A->rk;
    rank += B->rk;
    M = pastix_imax(M2, M1);
    N = pastix_imax(N2, N1);

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

    /*
     * A is rank null, nothing to do
     */
    if (A->rk == 0) {
        return rank;
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
        if ( A->rk == -1 ) {
            /*
             * TODO: This case can be improved by compressing A, and then
             * copying it into B, however the criteria to keep A compressed or
             * not must be based on B dimension, and not on A ones
             */
            u = malloc( M * N * sizeof(pastix_complex64_t) );

            if ( M1 != M || N1 != N ) {
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, N,
                                     0., 0., u, M );
            }
            ret = core_zgeadd( PastixNoTrans, M1, N1,
                               alpha, A->u, ldau,
                               0., u + M * offy + offx, M );
            assert(ret == 0);

            core_zge2lr_RRQR( tol, M, N, u, M, B );
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
    /* tauV */
    wzsize += minV;

    /* RRQR workspaces */
    ldwork = pastix_imax(rank, N);
    wzsize += (2 * nb + 1) * ldwork;

    zbuf = malloc( wzsize * sizeof(pastix_complex64_t) + 2 * pastix_imax(rank, N) * sizeof(double) );

    u1u2  = zbuf;
    v1v2  = u1u2 + M * rank;
    tauV  = v1v2 + N * rank;
    zwork = tauV + rank;

    rwork = (double*)(zbuf + wzsize);

    MALLOC_INTERN( jpvt, pastix_imax(rank, N), pastix_int_t );

    /*
     * Concatenate U2 and U1 in u1u2
     *  [ u2  0  ]
     *  [ u2  u1 ]
     *  [ u2  0  ]
     */
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, B->rk,
                         B->u, ldbu, u1u2, M );

    tmp = u1u2 + B->rk * M;
    if ( A->rk == -1 ) {
        /*
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
            /*
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
    /*
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

    /*
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
        /*
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
        /*
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
    /*
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


    /*
     * Perform RRQR factorization on v1v2 = (Q2 R2)
     */

    /*
     * Orthogonalize [u2, u1]
     * [u2, u1] = [u2, u1 - u2(u2Tu1)] * (I u2Tu1)
     *                                   (0   I  )
     */
    {
        pastix_int_t rA = A->rk;
        pastix_int_t rB = B->rk;
        if (A->rk == -1 && M1 >= N1)
            rA = N1;
        if (A->rk == -1 && M1 < N1)
            rA = M1;

        /* We do not care is A was integrated into v1v2 */
        if (rA != 0){
            pastix_complex64_t *u2Tu1 = malloc(rA * rB * sizeof(pastix_complex64_t));
            pastix_complex64_t *tmpU = u1u2 + offx;
            pastix_complex64_t *tmpV = v1v2 + rB;

            /* Form u2Tu1 */
            if (rA == N1){
                cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                            rB, rA, N1,
                            CBLAS_SADDR(zone), tmpU, ldbu,
                            A->u, ldau,
                            CBLAS_SADDR(zzero), u2Tu1, rB );
            }
            else if (rA == M1){
                pastix_int_t i, j;
                for (i=0; i<rA; i++){
                    for (j=0; j<rB; j++){
                        u2Tu1[rB * i + j] = tmpU[ldbu * j + i];
                    }
                }
            }
            else {
                cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                            rB, rA, rA,
                            CBLAS_SADDR(zone), tmpU, ldbu,
                            A->u, ldau,
                            CBLAS_SADDR(zzero), u2Tu1, rB );
            }

            /* Orthogonalize u1u2 */
            tmpU = u1u2 + B->rk * M;
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        M, rA, rB,
                        CBLAS_SADDR(mzone), u1u2, M,
                        u2Tu1, rB,
                        CBLAS_SADDR(zone), tmpU, M );

            /* Update v1v2 */
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        rB, N, rA,
                        CBLAS_SADDR(zone), u2Tu1, rB,
                        tmpV, rank,
                        CBLAS_SADDR(zone), v1v2, rank );

            /* Orthonormalize u1u2 */
            {
                pastix_int_t i, j;
                for (i=0; i<rank; i++){
                    pastix_complex64_t *tmpU = u1u2 + M * i;
                    pastix_complex64_t *tmpV = v1v2 + i;
                    double norm = cblas_dznrm2(M, tmpU, 1);

                    if (norm > LAPACKE_dlamch('e')){
                        for (j=0; j<M; j++){
                            tmpU[j] /= (norm);
                        }
                        for (j=0; j<N; j++){
                            tmpV[rank * j] *= norm;
                        }
                    }
                    else{
                        for (j=0; j<M; j++){
                            tmpU[j] = 0.;
                        }
                        for (j=0; j<N; j++){
                            tmpV[rank * j] = 0.;
                        }
                    }
                }
            }

            memFree_null(u2Tu1);
        }
    }

    norm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', rank, N,
                                v1v2, rank, NULL );

    new_rank = core_zrrqr(rank, N,
                          v1v2, rank,
                          jpvt, tauV,
                          zwork, ldwork,
                          rwork,
                          tol * norm, nb, rank-1);

    /*
     * First case: The rank is too big, so we decide to uncompress the result
     */
    if ( new_rank*2 > pastix_imin( M, N ) || new_rank == -1) {
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
        memFree_null(jpvt);
        return 0;
    }
    else if ( new_rank == 0 ) {
        core_zlrfree(B);
        memFree_null(zbuf);
        memFree_null(jpvt);
        return 0;
    }

    /*
     * We need to reallocate the buffer to store the new compressed version of B
     * because it wasn't big enough
     */
    ret = core_zlrsze( 0, M, N, B, new_rank, -1 );
    assert( ret != -1 );
    assert( B->rkmax >= new_rank );
    assert( B->rkmax >= B->rk    );

    ldbv = B->rkmax;

    pastix_complex64_t *work, *work2;
    work = malloc(N * new_rank * sizeof(pastix_complex64_t));
    work2 = malloc(N * new_rank * sizeof(pastix_complex64_t));

    memset(work, 0, N * new_rank * sizeof(pastix_complex64_t));

    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'U', new_rank, N,
                               v1v2, rank,
                               work, new_rank );
    assert(ret == 0);

    for (i=0; i<N; i++){
        memcpy(work2 + jpvt[i] * new_rank,
               work  + i * new_rank,
               new_rank * sizeof(pastix_complex64_t));
    }

    /* Compute Q2 factor */
    ret = LAPACKE_zungqr( LAPACK_COL_MAJOR, rank, new_rank, new_rank,
                          v1v2, rank, tauV );
    assert(ret == 0);

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                M, new_rank, rank,
                CBLAS_SADDR(zone),  u1u2, M,
                v1v2, rank,
                CBLAS_SADDR(zzero), B->u, ldbu);

    memcpy(B->v, work2, new_rank * N * sizeof(pastix_complex64_t));

    memFree_null(zbuf);
    memFree_null(jpvt);
    memFree_null(work);
    memFree_null(work2);

    (void)ret;
    return new_rank;
}

/**
 *******************************************************************************
 *
 * @brief Arithmetic free interface to core_zge2lr_RRQR
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
void
core_zge2lr_RRQR_interface( pastix_fixdbl_t tol, pastix_int_t m, pastix_int_t n,
                            const void *Aptr, pastix_int_t lda,
                            void *Alr )
{
    const pastix_complex64_t *A = (const pastix_complex64_t *) Aptr;
    core_zge2lr_RRQR( tol, m, n, A, lda, Alr );
}

/**
 *******************************************************************************
 *
 * @brief Arithmetic free interface to core_zrradd_RRQR
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          The absolute tolerance criteria
 *
 * @param[in] transA1
 *         @arg PastixNoTrans:  No transpose, op( A ) = A;
 *         @arg PastixTrans:  Transpose, op( A ) = A';
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
 * @return  The new rank of u2 v2^T or -1 if ranks are too large for
 *          recompression
 *
 *******************************************************************************/
int
core_zrradd_RRQR_interface( pastix_fixdbl_t tol, pastix_trans_t transA1, const void *alphaptr,
                            pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                            pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                            pastix_int_t offx, pastix_int_t offy )
{
    const pastix_complex64_t *alpha = (const pastix_complex64_t *) alphaptr;
    return core_zrradd_RRQR( tol, transA1, *alpha,
                             M1, N1, A,
                             M2, N2, B,
                             offx, offy );
}
