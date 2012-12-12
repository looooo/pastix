/**
 *
 * @file pzgemm.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m, n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define B(m, n) BLKADDR(B, PLASMA_Complex64_t, m, n)
#define C(m, n) BLKADDR(C, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel tile matrix-matrix multiplication - static scheduling
 **/
void plasma_pzgemm(plasma_context_t *plasma)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    PLASMA_Complex64_t alpha;
    PLASMA_desc A;
    PLASMA_desc B;
    PLASMA_Complex64_t beta;
    PLASMA_desc C;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int K, X, Y;
    int k, m, n;
    int next_m;
    int next_n;
    int ldam, ldak, ldbn, ldbk, ldcm;

    PLASMA_Complex64_t zbeta;
    PLASMA_Complex64_t zone = (PLASMA_Complex64_t)1.0;

    plasma_unpack_args_9(transA, transB, alpha, A, B, beta, C, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= C.mt && n < C.nt) {
        n++;
        m = m-C.mt;
    }

    while (n < C.nt) {
        next_m = m;
        next_n = n;

        next_m += PLASMA_SIZE;
        while (next_m >= C.mt && next_n < C.nt) {
            next_n++;
            next_m = next_m - C.mt;
        }

        X = m == C.mt-1 ? C.m - m*C.mb : C.mb;
        Y = n == C.nt-1 ? C.n - n*C.nb : C.nb;

        ldcm = BLKLDD(C, m);
        /*
         *  A: PlasmaNoTrans / B: PlasmaNoTrans
         */
        if (transA == PlasmaNoTrans) {
            ldam = BLKLDD(A, m);
            if (transB == PlasmaNoTrans) {
                for (k = 0; k < A.nt; k++) {
                    K = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                    ldbk = BLKLDD(B, k);
                    zbeta = k == 0 ? beta : zone;
                    CORE_zgemm(
                        transA, transB,
                        X, Y, K,
                        alpha, A(m, k), ldam,
                               B(k, n), ldbk,
                        zbeta, C(m, n), ldcm);
                }
            }
            /*
             *  A: PlasmaNoTrans / B: Plasma[Conj]Trans
             */
            else {
                ldbn = BLKLDD(B, n);
                for (k = 0; k < A.nt; k++) {
                    K = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                    zbeta = k == 0 ? beta : zone;
                    CORE_zgemm(
                        transA, transB,
                        X, Y, K,
                        alpha, A(m, k), ldam,
                               B(n, k), ldbn,
                        zbeta, C(m, n), ldcm);
                }
            }
        }
        /*
         *  A: Plasma[Conj]Trans / B: PlasmaNoTrans
         */
        else {
            if (transB == PlasmaNoTrans) {
                for (k = 0; k < A.mt; k++) {
                    K = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    zbeta = k == 0 ? beta : zone;
                    CORE_zgemm(
                        transA, transB,
                        X, Y, K,
                        alpha, A(k, m), ldak,
                               B(k, n), ldbk,
                        zbeta, C(m, n), ldcm);
                }
            }
            /*
             *  A: Plasma[Conj]Trans / B: Plasma[Conj]Trans
             */
            else {
                ldbn = BLKLDD(B, n);
                for (k = 0; k < A.mt; k++) {
                    K = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                    ldak = BLKLDD(A, k);
                    zbeta = k == 0 ? beta : zone;
                    CORE_zgemm(
                        transA, transB,
                        X, Y, K,
                        alpha, A(k, m), ldak,
                               B(n, k), ldbn,
                        zbeta, C(m, n), ldcm);
                }
            }
        }
        m = next_m;
        n = next_n;
    }
}

/***************************************************************************//**
 *  Parallel tile matrix-matrix multiplication - dynamic scheduling
 **/
void plasma_pzgemm_quark(PLASMA_enum transA, PLASMA_enum transB,
                         PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B,
                         PLASMA_Complex64_t beta,  PLASMA_desc C,
                         PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int m, n, k;
    int ldam, ldak, ldbn, ldbk, ldcm;
    int tempmm, tempnn, tempkn, tempkm;

    PLASMA_Complex64_t zbeta;
    PLASMA_Complex64_t zone = (PLASMA_Complex64_t)1.0;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (m = 0; m < C.mt; m++) {
        tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
        ldcm = BLKLDD(C, m);
        for (n = 0; n < C.nt; n++) {
            tempnn = n == C.nt-1 ? C.n-n*C.nb : C.nb;
            /*
             *  A: PlasmaNoTrans / B: PlasmaNoTrans
             */
            if (transA == PlasmaNoTrans) {
                ldam = BLKLDD(A, m);
                if (transB == PlasmaNoTrans) {
                    for (k = 0; k < A.nt; k++) {
                        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        QUARK_CORE_zgemm(
                            plasma->quark, &task_flags,
                            transA, transB,
                            tempmm, tempnn, tempkn, A.mb,
                            alpha, A(m, k), ldam,  /* lda * Z */
                                   B(k, n), ldbk,  /* ldb * Y */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
                /*
                 *  A: PlasmaNoTrans / B: Plasma[Conj]Trans
                 */
                else {
                    ldbn = BLKLDD(B, n);
                    for (k = 0; k < A.nt; k++) {
                        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                        zbeta = k == 0 ? beta : zone;
                        QUARK_CORE_zgemm(
                            plasma->quark, &task_flags,
                            transA, transB,
                            tempmm, tempnn, tempkn, A.mb,
                            alpha, A(m, k), ldam,  /* lda * Z */
                                   B(n, k), ldbn,  /* ldb * Z */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
            }
            /*
             *  A: Plasma[Conj]Trans / B: PlasmaNoTrans
             */
            else {
                if (transB == PlasmaNoTrans) {
                    for (k = 0; k < A.mt; k++) {
                        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        QUARK_CORE_zgemm(
                            plasma->quark, &task_flags,
                            transA, transB,
                            tempmm, tempnn, tempkm, A.mb,
                            alpha, A(k, m), ldak,  /* lda * X */
                                   B(k, n), ldbk,  /* ldb * Y */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
                /*
                 *  A: Plasma[Conj]Trans / B: Plasma[Conj]Trans
                 */
                else {
                    ldbn = BLKLDD(B, n);
                    for (k = 0; k < A.mt; k++) {
                        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        QUARK_CORE_zgemm(
                            plasma->quark, &task_flags,
                            transA, transB,
                            tempmm, tempnn, tempkm, A.mb,
                            alpha, A(k, m), ldak,  /* lda * X */
                                   B(n, k), ldbn,  /* ldb * Z */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
            }
        }
    }
}
