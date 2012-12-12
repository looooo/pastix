/**
 *
 * @file core_ztstrf.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <cblas.h>
#include <math.h>

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_ztstrf computes an LU factorization of a complex matrix formed
 *  by an upper triangular NB-by-N tile U on top of a M-by-N tile A
 *  using partial pivoting with row interchanges.
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] NB
 *
 * @param[in,out] U
 *         On entry, the NB-by-N upper triangular tile.
 *         On exit, the new factor U from the factorization
 *
 * @param[in] LDU
 *         The leading dimension of the array U.  LDU >= max(1,NB).
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile to be factored.
 *         On exit, the factor L from the factorization
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[in,out] L
 *         On entry, the IB-by-N lower triangular tile.
 *         On exit, the interchanged rows form the tile A in case of pivoting.
 *
 * @param[in] LDL
 *         The leading dimension of the array L.  LDL >= max(1,IB).
 *
 * @param[out] IPIV
 *         The pivot indices; for 1 <= i <= min(M,N), row i of the
 *         tile U was interchanged with row IPIV(i) of the tile A.
 *
 * @param[in,out] WORK
 *
 * @param[in] LDWORK
 *         The leading dimension of the array WORK.
 *
 * @param[out] INFO
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *         \retval >0 if INFO = k, U(k,k) is exactly zero. The factorization
 *              has been completed, but the factor U is exactly
 *              singular, and division by zero will occur if it is used
 *              to solve a system of equations.
 *
 ******************************************************************************/

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztstrf = PCORE_ztstrf
#define CORE_ztstrf PCORE_ztstrf
#define CORE_zssssm PCORE_zssssm
int  CORE_zssssm(int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex64_t *A1, int LDA1,
                 PLASMA_Complex64_t *A2, int LDA2,
                 const PLASMA_Complex64_t *L1, int LDL1,
                 const PLASMA_Complex64_t *L2, int LDL2,
                 const int *IPIV);
#endif
int CORE_ztstrf(int M, int N, int IB, int NB,
                PLASMA_Complex64_t *U, int LDU,
                PLASMA_Complex64_t *A, int LDA,
                PLASMA_Complex64_t *L, int LDL,
                int *IPIV,
                PLASMA_Complex64_t *WORK, int LDWORK,
                int *INFO)
{
    static PLASMA_Complex64_t zzero = 0.0;
    static PLASMA_Complex64_t mzone =-1.0;

    PLASMA_Complex64_t alpha;
    int i, j, ii, sb;
    int im, ip;

    /* Check input arguments */
    *INFO = 0;
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if (IB < 0) {
        coreblas_error(3, "Illegal value of IB");
        return -3;
    }
    if ((LDU < max(1,NB)) && (NB > 0)) {
        coreblas_error(6, "Illegal value of LDU");
        return -6;
    }
    if ((LDA < max(1,M)) && (M > 0)) {
        coreblas_error(8, "Illegal value of LDA");
        return -8;
    }
    if ((LDL < max(1,IB)) && (IB > 0)) {
        coreblas_error(10, "Illegal value of LDL");
        return -10;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    /* Set L to 0 */
    memset(L, 0, LDL*N*sizeof(PLASMA_Complex64_t));

    ip = 0;
    for (ii = 0; ii < N; ii += IB) {
        sb = min(N-ii, IB);

        for (i = 0; i < sb; i++) {
            im = cblas_izamax(M, &A[LDA*(ii+i)], 1);
            IPIV[ip] = ii+i+1;

            if (cabs(A[LDA*(ii+i)+im]) > cabs(U[LDU*(ii+i)+ii+i])) {
                /*
                 * Swap behind.
                 */
                cblas_zswap(i, &L[LDL*ii+i], LDL, &WORK[im], LDWORK );
                /*
                 * Swap ahead.
                 */
                cblas_zswap(sb-i, &U[LDU*(ii+i)+ii+i], LDU, &A[LDA*(ii+i)+im], LDA );
                /*
                 * Set IPIV.
                 */
                IPIV[ip] = NB + im + 1;

                for (j = 0; j < i; j++) {
                    A[LDA*(ii+j)+im] = zzero;
                }
            }

            if ((*INFO == 0) && (cabs(A[LDA*(ii+i)+im]) == zzero)
                && (cabs(U[LDU*(ii+i)+ii+i]) == zzero)) {
                *INFO = ii+i+1;
            }

            alpha = ((PLASMA_Complex64_t)1. / U[LDU*(ii+i)+ii+i]);
            cblas_zscal(M, CBLAS_SADDR(alpha), &A[LDA*(ii+i)], 1);
            cblas_zcopy(M, &A[LDA*(ii+i)], 1, &WORK[LDWORK*i], 1);
            cblas_zgeru(
                CblasColMajor, M, sb-i-1,
                CBLAS_SADDR(mzone), &A[LDA*(ii+i)], 1,
                &U[LDU*(ii+i+1)+ii+i], LDU,
                &A[LDA*(ii+i+1)], LDA );
            ip = ip+1;
        }
        /*
         * Apply the subpanel to the rest of the panel.
         */
        if(ii+i < N) {
            for(j = ii; j < ii+sb; j++) {
                if (IPIV[j] <= NB) {
                    IPIV[j] = IPIV[j] - ii;
                }
            }

            CORE_zssssm(
                NB, N-(ii+sb), M, N-(ii+sb), sb, sb,
                &U[LDU*(ii+sb)+ii], LDU,
                &A[LDA*(ii+sb)], LDA,
                &L[LDL*ii], LDL,
                WORK, LDWORK, &IPIV[ii]);

            for(j = ii; j < ii+sb; j++) {
                if (IPIV[j] <= NB) {
                    IPIV[j] = IPIV[j] + ii;
                }
            }
        }
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ztstrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex64_t *U, int ldu,
                       PLASMA_Complex64_t *A, int lda,
                       PLASMA_Complex64_t *L, int ldl,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo)
{
    DAG_CORE_TSTRF;
    QUARK_Insert_Task(quark, CORE_ztstrf_quark, task_flags,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(int),                        &nb,            VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    U,                     INOUT | QUARK_REGION_D | QUARK_REGION_U,
        sizeof(int),                        &ldu,           VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,                     INOUT | LOCALITY,
        sizeof(int),                        &lda,           VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    L,                     OUTPUT,
        sizeof(int),                        &ldl,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(PLASMA_Complex64_t)*ib*nb,    NULL,                  SCRATCH,
        sizeof(int),                        &nb,            VALUE,
        sizeof(PLASMA_sequence*),           &sequence,      VALUE,
        sizeof(PLASMA_request*),            &request,       VALUE,
        sizeof(PLASMA_bool),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztstrf_quark = PCORE_ztstrf_quark
#define CORE_ztstrf_quark PCORE_ztstrf_quark
#endif
void CORE_ztstrf_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    int nb;
    PLASMA_Complex64_t *U;
    int ldu;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_Complex64_t *L;
    int ldl;
    int *IPIV;
    PLASMA_Complex64_t *WORK;
    int ldwork;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_bool check_info;
    int iinfo;

    int info;

    quark_unpack_args_17(quark, m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK, ldwork, sequence, request, check_info, iinfo);
    CORE_ztstrf(m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK, ldwork, &info);
    if (info != PLASMA_SUCCESS && check_info)
        plasma_sequence_flush(quark, sequence, request, iinfo + info);
}
