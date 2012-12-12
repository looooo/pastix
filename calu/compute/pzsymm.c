/**
 *
 * @file pzsymm.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define B(m,n) BLKADDR(B, PLASMA_Complex64_t, m, n)
#define C(m,n) BLKADDR(C, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel tile symmetric matrix-matrix multiplication - static scheduling
 **/
void plasma_pzsymm(plasma_context_t *plasma)
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
    int lda, ldak, ldb, ldc;
    int tempmm, tempnn, tempkn, tempkm;

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

        ldc = BLKLDD(C, m);
        /*
         *  PlasmaLeft / PlasmaLower
         */
        if (side == PlasmaLeft) {
            lda = BLKLDD(A, m);
            if (uplo == PlasmaLower) {
                for (k = 0; k < C.mt; k++) {
                    tempkm = k == C.mt-1 ? C.m-k*C.mb : C.mb;
                    ldak = BLKLDD(A, k);
                    ldb  = BLKLDD(B, k);
                    zbeta = k == 0 ? beta : zone;
                    if (k < m) {
                        CORE_zgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm,
                            alpha, A(m, k), lda,
                                   B(k, n), ldb,
                            zbeta, C(m, n), ldc);
                    }
                    else {
                        if (k == m) {
                            CORE_zsymm(
                                side, uplo,
                                tempmm, tempnn,
                                alpha, A(k, k), ldak,
                                       B(k, n), ldb,
                                zbeta, C(m, n), ldc);
                        }
                        else {
                            CORE_zgemm(
                                PlasmaTrans, PlasmaNoTrans,
                                tempmm, tempnn, tempkm,
                                alpha, A(k, m), ldak,
                                       B(k, n), ldb,
                                zbeta, C(m, n), ldc);
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
                    ldb  = BLKLDD(B, k);
                    zbeta = k == 0 ? beta : zone;
                    if (k < m) {
                        CORE_zgemm(
                            PlasmaTrans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm,
                            alpha, A(k, m), ldak,
                                   B(k, n), ldb,
                            zbeta, C(m, n), ldc);
                    }
                    else {
                        if (k == m) {
                            CORE_zsymm(
                                side, uplo,
                                tempmm, tempnn,
                                alpha, A(k, k), ldak,
                                       B(k, n), ldb,
                                zbeta, C(m, n), ldc);
                        }
                        else {
                            CORE_zgemm(
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, tempkm,
                                alpha, A(m, k), lda,
                                       B(k, n), ldb,
                                zbeta, C(m, n), ldc);
                        }
                    }
                }
            }
        }
        /*
         *  PlasmaRight / PlasmaLower
         */
        else {
            lda = BLKLDD(A, n);
            ldb = BLKLDD(B, m);
            if (uplo == PlasmaLower) {
                for (k = 0; k < C.nt; k++) {
                    tempkn = k == C.nt-1 ? C.n-k*C.nb : C.nb;
                    ldak = BLKLDD(A, k);
                    zbeta = k == 0 ? beta : zone;
                    if (k < n) {
                        CORE_zgemm(
                            PlasmaNoTrans, PlasmaTrans,
                            tempmm, tempnn, tempkn,
                            alpha, B(m, k), ldb,
                                   A(n, k), lda,
                            zbeta, C(m, n), ldc);
                    }
                    else {
                        if (n == k) {
                            CORE_zsymm(
                                side, uplo,
                                tempmm, tempnn,
                                alpha, A(k, k), ldak,
                                       B(m, k), ldb,
                                zbeta, C(m, n), ldc);
                        }
                        else {
                            CORE_zgemm(
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, tempkn,
                                alpha, B(m, k), ldb,
                                       A(k, n), ldak,
                                zbeta, C(m, n), ldc);
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
                            alpha, B(m, k), ldb,
                                   A(k, n), ldak,
                            zbeta, C(m, n), ldc);
                    }
                    else {
                        if (n == k) {
                            CORE_zsymm(
                                side, uplo,
                                tempmm, tempnn,
                                alpha, A(k, k), ldak,
                                       B(m, k), ldb,
                                zbeta, C(m, n), ldc);
                        }
                        else {
                            CORE_zgemm(
                                PlasmaNoTrans, PlasmaTrans,
                                tempmm, tempnn, tempkn,
                                alpha, B(m, k), ldb,
                                       A(n, k), lda,
                                zbeta, C(m, n), ldc);
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
 *  Parallel tile symmetric matrix-matrix multiplication - dynamic scheduling
 **/
void plasma_pzsymm_quark(PLASMA_enum side, PLASMA_enum uplo,
                          PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B,
                          PLASMA_Complex64_t beta, PLASMA_desc C,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int lda, ldak, ldb, ldc;
    int tempmm, tempnn, tempkn, tempkm;

    PLASMA_Complex64_t zbeta;
    PLASMA_Complex64_t zone = (PLASMA_Complex64_t)1.0;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (m = 0; m < C.mt; m++) {
        tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
        ldc = BLKLDD(C, m);
        for (n = 0; n < C.nt; n++) {
            tempnn = n == C.nt-1 ? C.n-n*C.nb : C.nb;
            /*
             *  PlasmaLeft / PlasmaLower
             */
            if (side == PlasmaLeft) {
                lda = BLKLDD(A, m);
                if (uplo == PlasmaLower) {
                    for (k = 0; k < C.mt; k++) {
                        tempkm = k == C.mt-1 ? C.m-k*C.mb : C.mb;
                        ldak = BLKLDD(A, k);
                        ldb  = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < m) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, tempkm, A.mb,
                                alpha, A(m, k), lda,  /* lda * K */
                                       B(k, n), ldb,  /* ldb * Y */
                                zbeta, C(m, n), ldc); /* ldc * Y */
                        }
                        else {
                            if (k == m) {
                                QUARK_CORE_zsymm(
                                    plasma->quark, &task_flags,
                                    side, uplo,
                                    tempmm, tempnn, A.mb,
                                    alpha, A(k, k), ldak, /* ldak * X */
                                           B(k, n), ldb,  /* ldb  * Y */
                                    zbeta, C(m, n), ldc); /* ldc  * Y */
                            }
                            else {
                                QUARK_CORE_zgemm(
                                    plasma->quark, &task_flags,
                                    PlasmaTrans, PlasmaNoTrans,
                                    tempmm, tempnn, tempkm, A.mb,
                                    alpha, A(k, m), ldak, /* ldak * X */
                                           B(k, n), ldb,  /* ldb  * Y */
                                    zbeta, C(m, n), ldc); /* ldc  * Y */
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
                        ldb  = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < m) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaTrans, PlasmaNoTrans,
                                tempmm, tempnn, tempkm, A.mb,
                                alpha, A(k, m), ldak, /* ldak * X */
                                       B(k, n), ldb,  /* ldb  * Y */
                                zbeta, C(m, n), ldc); /* ldc  * Y */
                        }
                        else {
                            if (k == m) {
                                QUARK_CORE_zsymm(
                                    plasma->quark, &task_flags,
                                    side, uplo,
                                    tempmm, tempnn, A.mb,
                                    alpha, A(k, k), ldak, /* ldak * K */
                                           B(k, n), ldb,  /* ldb  * Y */
                                    zbeta, C(m, n), ldc); /* ldc  * Y */
                            }
                            else {
                                QUARK_CORE_zgemm(
                                    plasma->quark, &task_flags,
                                    PlasmaNoTrans, PlasmaNoTrans,
                                    tempmm, tempnn, tempkm, A.mb,
                                    alpha, A(m, k), lda,  /* lda * K */
                                           B(k, n), ldb,  /* ldb * Y */
                                    zbeta, C(m, n), ldc); /* ldc * Y */
                            }
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaLower
             */
            else {
                lda = BLKLDD(A, n);
                ldb = BLKLDD(B, m);
                if (uplo == PlasmaLower) {
                    for (k = 0; k < C.nt; k++) {
                        tempkn = k == C.nt-1 ? C.n-k*C.nb : C.nb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < n) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaTrans,
                                tempmm, tempnn, tempkn, A.mb,
                                alpha, B(m, k), ldb,  /* ldb * K */
                                       A(n, k), lda,  /* lda * K */
                                zbeta, C(m, n), ldc); /* ldc * Y */
                        }
                        else {
                            if (k == n) {
                                QUARK_CORE_zsymm(
                                    plasma->quark, &task_flags,
                                    side, uplo,
                                    tempmm, tempnn, A.mb,
                                    alpha, A(k, k), ldak, /* ldak * Y */
                                           B(m, k), ldb,  /* ldb  * Y */
                                    zbeta, C(m, n), ldc); /* ldc  * Y */
                            }
                            else {
                                QUARK_CORE_zgemm(
                                    plasma->quark, &task_flags,
                                    PlasmaNoTrans, PlasmaNoTrans,
                                    tempmm, tempnn, tempkn, A.mb,
                                    alpha, B(m, k), ldb,  /* ldb  * K */
                                           A(k, n), ldak, /* ldak * Y */
                                    zbeta, C(m, n), ldc); /* ldc  * Y */
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
                                alpha, B(m, k), ldb,  /* ldb  * K */
                                       A(k, n), ldak, /* ldak * Y */
                                zbeta, C(m, n), ldc); /* ldc  * Y */
                        }
                        else {
                            if (k == n) {
                                QUARK_CORE_zsymm(
                                    plasma->quark, &task_flags,
                                    side, uplo,
                                    tempmm, tempnn, A.mb,
                                    alpha, A(k, k), ldak, /* ldak * Y */
                                           B(m, k), ldb,  /* ldb  * Y */
                                    zbeta, C(m, n), ldc); /* ldc  * Y */
                            }
                            else {
                                QUARK_CORE_zgemm(
                                    plasma->quark, &task_flags,
                                    PlasmaNoTrans, PlasmaTrans,
                                    tempmm, tempnn, tempkn, A.mb,
                                    alpha, B(m, k), ldb,  /* ldb * K */
                                           A(n, k), lda,  /* lda * K */
                                    zbeta, C(m, n), ldc); /* ldc * Y */
                            }
                        }
                    }
                }
            }
        }
    }
}
