/**
 *
 * @file pzherk.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define C(m,n) BLKADDR(C, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel tile Hermitian rank-k update - static scheduling
 **/
void plasma_pzherk(plasma_context_t *plasma)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    double alpha;
    PLASMA_desc A;
    double beta;
    PLASMA_desc C;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int m, n, k;
    int next_m;
    int next_n;
    int ldam, ldan, ldak, ldcm, ldcn;
    int tempkn, tempkm, tempmm, tempnn;

    PLASMA_Complex64_t zone = (PLASMA_Complex64_t)1.0;
    PLASMA_Complex64_t zalpha;
    PLASMA_Complex64_t zbeta;
    double dbeta;

    plasma_unpack_args_8(uplo, trans, alpha, A, beta, C, sequence, request);
    zalpha = (PLASMA_Complex64_t)alpha;
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= C.mt && n < C.nt) {
        n++;
        m = m-C.mt+n;
    }

    while (n < C.nt) {
        next_n = n;
        next_m = m + PLASMA_SIZE;
        while (next_m >= C.mt && next_n < C.nt) {
            next_n++;
            next_m = next_m - C.mt + next_n;
        }

        tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
        tempnn = n == C.nt-1 ? C.n-n*C.nb : C.nb;

        if (m == n) {
            ldcm = BLKLDD(C, m);
            /*
             *  PlasmaNoTrans
             */
            if (trans == PlasmaNoTrans) {
                ldam = BLKLDD(A, m);
                for (k = 0; k < A.nt; k++) {
                    tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                    dbeta = k == 0 ? beta : 1.0;
                    CORE_zherk(
                        uplo, trans,
                        tempnn, tempkn,
                        alpha, A(m, k), ldam,
                        dbeta, C(m, n), ldcm);
                }
            }
            /*
             *  Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < A.mt; k++) {
                    tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                    ldak = BLKLDD(A, k);
                    dbeta = k == 0 ? beta : 1.0;
                    CORE_zherk(
                        uplo, trans,
                        tempnn, tempkm,
                        alpha, A(k, m), ldak,
                        dbeta, C(m, n), ldcm);
                }
            }
        }
        else {
            if (trans == PlasmaNoTrans) {
                ldam = BLKLDD(A, m);
                ldan = BLKLDD(A, n);
                /*
                 *  PlasmaNoTrans / PlasmaLower
                 */
                if (uplo == PlasmaLower) {
                    ldcm = BLKLDD(C, m);
                    for (k = 0; k < A.nt; k++) {
                        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                        zbeta = k == 0 ? (PLASMA_Complex64_t)beta : zone;
                        CORE_zgemm(
                            trans, PlasmaConjTrans,
                            tempmm, tempnn, tempkn,
                            zalpha, A(m, k), ldam,
                                    A(n, k), ldan,
                            zbeta,  C(m, n), ldcm);
                    }
                }
                /*
                 *  PlasmaNoTrans / PlasmaUpper
                 */
                else {
                    ldcn = BLKLDD(C, n);
                    for (k = 0; k < A.nt; k++) {
                        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                        zbeta = k == 0 ? (PLASMA_Complex64_t)beta : zone;
                        CORE_zgemm(
                            trans, PlasmaConjTrans,
                            tempnn, tempmm, tempkn,
                            zalpha, A(n, k), ldan,
                                    A(m, k), ldam,
                            zbeta,  C(n, m), ldcn);
                    }
                }
            }
            else {
                /*
                 *  Plasma[Conj]Trans / PlasmaLower
                 */
                if (uplo == PlasmaLower) {
                    ldcm = BLKLDD(C, m);
                    for (k = 0; k < A.mt; k++) {
                        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? (PLASMA_Complex64_t)beta : zone;
                        CORE_zgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm,
                            zalpha, A(k, m), ldak,
                                    A(k, n), ldak,
                            zbeta,  C(m, n), ldcm);
                    }
                }
                /*
                 *  Plasma[Conj]Trans / PlasmaUpper
                 */
                else {
                    ldcn = BLKLDD(C, n);
                    for (k = 0; k < A.mt; k++) {
                        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? (PLASMA_Complex64_t)beta : zone;
                        CORE_zgemm(
                            trans, PlasmaNoTrans,
                            tempnn, tempmm, tempkm,
                            zalpha, A(k, n), ldak,
                                    A(k, m), ldak,
                            zbeta,  C(n, m), ldcn);
                    }
                }
            }
        }
        m = next_m;
        n = next_n;
    }
}

/***************************************************************************//**
 *  Parallel tile Hermitian rank-k update - dynamic scheduling
 **/
void plasma_pzherk_quark(PLASMA_enum uplo, PLASMA_enum trans,
                         double alpha, PLASMA_desc A,
                         double beta,  PLASMA_desc C,
                         PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int m, n, k;
    int ldak, ldam, ldan, ldcm, ldcn;
    int tempnn, tempmm, tempkn, tempkm;

    PLASMA_Complex64_t zone   = (PLASMA_Complex64_t)1.0;
    PLASMA_Complex64_t zalpha = (PLASMA_Complex64_t)alpha;
    PLASMA_Complex64_t zbeta;
    double dbeta;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (n = 0; n < C.nt; n++) {
        tempnn = n == C.nt-1 ? C.n-n*C.nb : C.nb;
        ldan = BLKLDD(A, n);
        ldcn = BLKLDD(C, n);
        /*
         *  PlasmaNoTrans
         */
        if (trans == PlasmaNoTrans) {
            for (k = 0; k < A.nt; k++) {
                tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                dbeta = k == 0 ? beta : 1.0;
                QUARK_CORE_zherk(
                    plasma->quark, &task_flags,
                    uplo, trans,
                    tempnn, tempkn, A.mb,
                    alpha, A(n, k), ldan, /* ldan * K */
                    dbeta, C(n, n), ldcn); /* ldc  * N */
            }
            /*
             *  PlasmaNoTrans / PlasmaLower
             */
            if (uplo == PlasmaLower) {
                for (m = n+1; m < C.mt; m++) {
                    tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
                    ldam = BLKLDD(A, m);
                    ldcm = BLKLDD(C, m);
                    for (k = 0; k < A.nt; k++) {
                        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                        zbeta = k == 0 ? (PLASMA_Complex64_t)beta : zone;
                        QUARK_CORE_zgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaConjTrans,
                            tempmm, tempnn, tempkn, A.mb,
                            zalpha, A(m, k), ldam,  /* ldam * K */
                                    A(n, k), ldan,  /* ldan * K */
                            zbeta,  C(m, n), ldcm); /* ldc  * N */
                    }
                }
            }
            /*
             *  PlasmaNoTrans / PlasmaUpper
             */
            else {
                for (m = n+1; m < C.mt; m++) {
                    tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
                    ldam = BLKLDD(A, m);
                    for (k = 0; k < A.nt; k++) {
                        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                        zbeta = k == 0 ? (PLASMA_Complex64_t)beta : zone;
                        QUARK_CORE_zgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaConjTrans,
                            tempnn, tempmm, tempkn, A.mb,
                            zalpha, A(n, k), ldan,  /* ldan * K */
                                    A(m, k), ldam,  /* ldam * M */
                            zbeta,  C(n, m), ldcn); /* ldc  * M */
                    }
                }
            }
        }
        /*
         *  Plasma[Conj]Trans
         */
        else {
            for (k = 0; k < A.mt; k++) {
                tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                ldak = BLKLDD(A, k);
                dbeta = k == 0 ? beta : 1.0;
                QUARK_CORE_zherk(
                    plasma->quark, &task_flags,
                    uplo, trans,
                    tempnn, tempkm, A.mb,
                    alpha, A(k, n), ldak,  /* lda * N */
                    dbeta, C(n, n), ldcn); /* ldc * N */
            }
            /*
             *  Plasma[Conj]Trans / PlasmaLower
             */
            if (uplo == PlasmaLower) {
                for (m = n+1; m < C.mt; m++) {
                    tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
                    ldcm = BLKLDD(C, m);
                    for (k = 0; k < A.mt; k++) {
                        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? (PLASMA_Complex64_t)beta : zone;
                        QUARK_CORE_zgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm, A.mb,
                            zalpha, A(k, m), ldak,  /* lda * M */
                                    A(k, n), ldak,  /* lda * N */
                            zbeta,  C(m, n), ldcm); /* ldc * N */
                    }
                }
            }
            /*
             *  Plasma[Conj]Trans / PlasmaUpper
             */
            else {
                for (m = n+1; m < C.mt; m++) {
                    tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
                    for (k = 0; k < A.mt; k++) {
                        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? (PLASMA_Complex64_t)beta : zone;
                        QUARK_CORE_zgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaNoTrans,
                            tempnn, tempmm, tempkm, A.mb,
                            zalpha, A(k, n), ldak,  /* lda * K */
                                    A(k, m), ldak,  /* lda * M */
                            zbeta,  C(n, m), ldcn); /* ldc * M */
                    }
                }
            }
        }
    }
}
