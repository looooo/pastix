/**
 *
 * @file pzhemm.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define B(m,n) BLKADDR(B, PLASMA_Complex64_t, m, n)
#define C(m,n) BLKADDR(C, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel tile Hermitian matrix-matrix multiplication - static scheduling
 **/
void plasma_pzhemm(plasma_context_t *plasma)
{
    PLASMA_enum side;
    PLASMA_enum uplo;
    PLASMA_Complex64_t alpha;
    PLASMA_desc A;
    PLASMA_desc B;
    PLASMA_Complex64_t beta;
    PLASMA_desc C;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k, m, n;
    int next_m;
    int next_n;
    int ldan, ldak, ldam, ldbk, ldbm, ldcm;
    int tempmm, tempnn, tempkm, tempkn;

    PLASMA_Complex64_t zbeta;
    PLASMA_Complex64_t zone = (PLASMA_Complex64_t)1.0;

    plasma_unpack_args_9(side, uplo, alpha, A, B, beta, C, sequence, request);
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

        tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
        tempnn = n == C.nt-1 ? C.n-n*C.nb : C.nb;
        ldcm = BLKLDD(C, m);
        /*
         *  PlasmaLeft / PlasmaLower
         */
        if (side == PlasmaLeft) {
            ldam = BLKLDD(A, m);
            if (uplo == PlasmaLower) {
                for (k = 0; k < C.mt; k++) {
                    tempkm = k == C.mt-1 ? C.m-k*C.mb : C.mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    zbeta = k == 0 ? beta : zone;
                    if (k < m) {
                        CORE_zgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm,
                            alpha, A(m, k), ldam,
                                   B(k, n), ldbk,
                            zbeta, C(m, n), ldcm);
                    }
                    else {
                        if (k == m) {
                            CORE_zhemm(
                                side, uplo,
                                tempmm, tempnn,
                                alpha, A(k, k), ldak,
                                       B(k, n), ldbk,
                                zbeta, C(m, n), ldcm);
                        }
                        else {
                            CORE_zgemm(
                                PlasmaConjTrans, PlasmaNoTrans,
                                tempmm, tempnn, tempkm,
                                alpha, A(k, m), ldak,
                                       B(k, n), ldbk,
                                zbeta, C(m, n), ldcm);
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaUpper
             */
            else {
                for (k = 0; k < C.mt; k++) {
                    tempkm = k == C.mt-1 ? C.m-k*C.mb : C.mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    zbeta = k == 0 ? beta : zone;
                    if (k < m) {
                        CORE_zgemm(
                            PlasmaConjTrans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm,
                            alpha, A(k, m), ldak,
                                   B(k, n), ldbk,
                            zbeta, C(m, n), ldcm);
                    }
                    else {
                        if (k == m) {
                            CORE_zhemm(
                                side, uplo,
                                tempmm, tempnn,
                                alpha, A(k, k), ldak,
                                       B(k, n), ldbk,
                                zbeta, C(m, n), ldcm);
                        }
                        else {
                            CORE_zgemm(
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, tempkm,
                                alpha, A(m, k), ldam,
                                       B(k, n), ldbk,
                                zbeta, C(m, n), ldcm);
                        }
                    }
                }
            }
        }
        /*
         *  PlasmaRight / PlasmaLower
         */
        else {
            ldan = BLKLDD(A, n);
            ldbm = BLKLDD(B, m);
            if (uplo == PlasmaLower) {
                for (k = 0; k < C.nt; k++) {
                    tempkn = k == C.nt-1 ? C.n-k*C.nb : C.nb;
                    ldak = BLKLDD(A, k);
                    zbeta = k == 0 ? beta : zone;
                    if (k < n) {
                        CORE_zgemm(
                            PlasmaNoTrans, PlasmaConjTrans,
                            tempmm, tempnn, tempkn,
                            alpha, B(m, k), ldbm,
                                   A(n, k), ldan,
                            zbeta, C(m, n), ldcm);
                    }
                    else {
                        if (n == k) {
                            CORE_zhemm(
                                side, uplo,
                                tempmm, tempnn,
                                alpha, A(k, k), ldak,
                                       B(m, k), ldbm,
                                zbeta, C(m, n), ldcm);
                        }
                        else {
                            CORE_zgemm(
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, tempkn,
                                alpha, B(m, k), ldbm,
                                       A(k, n), ldak,
                                zbeta, C(m, n), ldcm);
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaUpper
             */
            else {
                for (k = 0; k < C.nt; k++) {
                    tempkn = k == C.nt-1 ? C.n-k*C.nb : C.nb;
                    ldak = BLKLDD(A, k);
                    zbeta = k == 0 ? beta : zone;
                    if (k < n) {
                        CORE_zgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            tempmm, tempnn, tempkn,
                            alpha, B(m, k), ldbm,
                                   A(k, n), ldak,
                            zbeta, C(m, n), ldcm);
                    }
                    else {
                        if (n == k) {
                            CORE_zhemm(
                                side, uplo,
                                tempmm, tempnn,
                                alpha, A(k, k), ldak,
                                       B(m, k), ldbm,
                                zbeta, C(m, n), ldcm);
                        }
                        else {
                            CORE_zgemm(
                                PlasmaNoTrans, PlasmaConjTrans,
                                tempmm, tempnn, tempkn,
                                alpha, B(m, k), ldbm,
                                       A(n, k), ldan,
                                zbeta, C(m, n), ldcm);
                        }
                    }
                }
            }
        }
        m = next_m;
        n = next_n;
    }
}

/***************************************************************************//**
 *  Parallel tile Hermitian matrix-matrix multiplication - dynamic scheduling
 **/
void plasma_pzhemm_quark(PLASMA_enum side, PLASMA_enum uplo,
                         PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B,
                         PLASMA_Complex64_t beta, PLASMA_desc C,
                         PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldam, ldan, ldak, ldbk, ldbm, ldcm;
    int tempmm, tempnn, tempkn, tempkm;

    PLASMA_Complex64_t zbeta;
    PLASMA_Complex64_t zone = (PLASMA_Complex64_t)1.0;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for(m = 0; m < C.mt; m++) {
        tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
        ldcm = BLKLDD(C, m);
        for(n = 0; n < C.nt; n++) {
            tempnn = n == C.nt-1 ? C.n-n*C.nb : C.nb;
            /*
             *  PlasmaLeft / PlasmaLower
             */
            if (side == PlasmaLeft) {
                ldam = BLKLDD(A, m);
                if (uplo == PlasmaLower) {
                    for (k = 0; k < C.mt; k++) {
                        tempkm = k == C.mt-1 ? C.m-k*C.mb : C.mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < m) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, tempkm, A.mb,
                                alpha, A(m, k), ldam,  /* lda * K */
                                       B(k, n), ldbk,  /* ldb * Y */
                                zbeta, C(m, n), ldcm); /* ldc * Y */
                        }
                        else {
                            if (k == m) {
                                QUARK_CORE_zhemm(
                                    plasma->quark, &task_flags,
                                    side, uplo,
                                    tempmm, tempnn, A.mb,
                                    alpha, A(k, k), ldak,  /* ldak * X */
                                           B(k, n), ldbk,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                QUARK_CORE_zgemm(
                                    plasma->quark, &task_flags,
                                    PlasmaConjTrans, PlasmaNoTrans,
                                    tempmm, tempnn, tempkm, A.mb,
                                    alpha, A(k, m), ldak,  /* ldak * X */
                                           B(k, n), ldbk,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                        }
                    }
                }
                /*
                 *  PlasmaLeft / PlasmaUpper
                 */
                else {
                    for (k = 0; k < C.mt; k++) {
                        tempkm = k == C.mt-1 ? C.m-k*C.mb : C.mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < m) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaConjTrans, PlasmaNoTrans,
                                tempmm, tempnn, tempkm, A.mb,
                                alpha, A(k, m), ldak,  /* ldak * X */
                                       B(k, n), ldbk,  /* ldb  * Y */
                                zbeta, C(m, n), ldcm); /* ldc  * Y */
                        }
                        else {
                            if (k == m) {
                                QUARK_CORE_zhemm(
                                    plasma->quark, &task_flags,
                                    side, uplo,
                                    tempmm, tempnn, A.mb,
                                    alpha, A(k, k), ldak,  /* ldak * K */
                                           B(k, n), ldbk,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                QUARK_CORE_zgemm(
                                    plasma->quark, &task_flags,
                                    PlasmaNoTrans, PlasmaNoTrans,
                                    tempmm, tempnn, tempkm, A.mb,
                                    alpha, A(m, k), ldam,  /* lda * K */
                                           B(k, n), ldbk,  /* ldb * Y */
                                    zbeta, C(m, n), ldcm); /* ldc * Y */
                            }
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaLower
             */
            else {
                ldan = BLKLDD(A, n);
                ldbm = BLKLDD(B, m);
                if (uplo == PlasmaLower) {
                    for (k = 0; k < C.nt; k++) {
                        tempkn = k == C.nt-1 ? C.n-k*C.nb : C.nb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < n) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaConjTrans,
                                tempmm, tempnn, tempkn, A.mb,
                                alpha, B(m, k), ldbm,  /* ldb * K */
                                       A(n, k), ldan,  /* lda * K */
                                zbeta, C(m, n), ldcm); /* ldc * Y */
                        }
                        else {
                            if (k == n) {
                                QUARK_CORE_zhemm(
                                    plasma->quark, &task_flags,
                                    side, uplo,
                                    tempmm, tempnn, A.mb,
                                    alpha, A(k, k), ldak,  /* ldak * Y */
                                           B(m, k), ldbm,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                QUARK_CORE_zgemm(
                                    plasma->quark, &task_flags,
                                    PlasmaNoTrans, PlasmaNoTrans,
                                    tempmm, tempnn, tempkn, A.mb,
                                    alpha, B(m, k), ldbm,  /* ldb  * K */
                                           A(k, n), ldak,  /* ldak * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                        }
                    }
                }
                /*
                 *  PlasmaRight / PlasmaUpper
                 */
                else {
                    for (k = 0; k < C.nt; k++) {
                        tempkn = k == C.nt-1 ? C.n-k*C.nb : C.nb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < n) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, tempkn, A.mb,
                                alpha, B(m, k), ldbm,  /* ldb  * K */
                                       A(k, n), ldak,  /* ldak * Y */
                                zbeta, C(m, n), ldcm); /* ldc  * Y */
                        }
                        else {
                            if (k == n) {
                                QUARK_CORE_zhemm(
                                    plasma->quark, &task_flags,
                                    side, uplo,
                                    tempmm, tempnn, A.mb,
                                    alpha, A(k, k), ldak,  /* ldak * Y */
                                           B(m, k), ldbm,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                QUARK_CORE_zgemm(
                                    plasma->quark, &task_flags,
                                    PlasmaNoTrans, PlasmaConjTrans,
                                    tempmm, tempnn, tempkn, A.mb,
                                    alpha, B(m, k), ldbm,  /* ldb * K */
                                           A(n, k), ldan,  /* lda * K */
                                    zbeta, C(m, n), ldcm); /* ldc * Y */
                            }
                        }
                    }
                }
            }
        }
    }
}
