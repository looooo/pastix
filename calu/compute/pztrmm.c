/**
 *
 * @file pztrmm.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define B(m,n) BLKADDR(B, PLASMA_Complex64_t, m, n)

#if 0
/***************************************************************************//**
 *  Parallel tile triangular matrix-matrix multiplication - static scheduling
 **/
void plasma_pztrmm(plasma_context_t *plasma)
{
    PLASMA_enum side;
    PLASMA_enum uplo;
    PLASMA_enum trans;
    PLASMA_enum diag;
    PLASMA_Complex64_t alpha;
    PLASMA_desc A;
    PLASMA_desc B;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k, m, n;
    int next_m;
    int next_n;
    int lda, ldak, ldb, ldbk;
    int tempkm, tempkn, tempmm, tempnn;

    PLASMA_Complex64_t *lB;
    PLASMA_Complex64_t zone  = (PLASMA_Complex64_t)1.0;

    plasma_unpack_args_9(side, uplo, trans, diag, alpha, A, B, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= B.mt && n < B.nt) {
        n++;
        m = m-B.mt;
    }

    while (n < B.nt) {
        next_m = m;
        next_n = n;

        next_m += PLASMA_SIZE;
        while (next_m >= B.mt && next_n < B.nt) {
            next_n++;
            next_m = next_m - B.mt;
        }

        tempmm = m == B.mt-1 ? B.m - m * B.mb : B.mb;
        tempnn = n == B.nt-1 ? B.n - n * B.nb : B.nb;
        ldb = m < B.lm1 ? B.mb : B.lm%B.mb;
        lB  = B(m, n);

        if ( side == PlasmaLeft ) {
            if ( uplo == PlasmaUpper ) {
                if ( trans == PlasmaNoTrans ) {
                    lda = m < A.lm1 ? A.mb : A.lm%B.mb;

                    CORE_ztrmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(m, m), lda,  /* lda * tempkm */
                               lB,      ldb); /* ldb * tempnn */

                    for (k = m+1; k < A.mt; k++) {
                        tempkn = k == A.nt-1 ? A.n - k * A.nb : A.nb;
                        ldbk = BLKLDD(B, k);

                        CORE_zgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, tempkn,
                            alpha, A(m, k), lda,
                                   B(k, n), ldbk,
                            zone,  lB,      ldb);
                    }
                }
                /*
                 *  PlasmaLeft / PlasmaUpper / Plasma[Conj]Trans
                 */
                else {
                    lda = m < A.lm1 ? A.mb : A.lm%B.mb;

                    CORE_ztrmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(m, m), lda,  /* lda * tempkm */
                               lB,      ldb); /* ldb * tempnn */

                    for (k = 0; k < m; k++) {
                        CORE_zgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, B.mb,
                            alpha, A(k, m), A.mb,
                                   B(k, n), B.mb,
                            zone,  lB,      ldb);
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaLower / PlasmaNoTrans
             */
            else {
                if ( trans == PlasmaNoTrans ) {
                    lda = m < A.lm1 ? A.mb : A.lm%B.mb;

                    CORE_ztrmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(m, m), lda,  /* lda * tempkm */
                               lB,      ldb); /* ldb * tempnn */

                    for (k = 0; k < m; k++) {
                        CORE_zgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, B.mb,
                            alpha, A(m, k), lda,
                                   B(k, n), B.mb,
                            zone,  lB,      ldb);
                    }
                }
                /*
                 *  PlasmaLeft / PlasmaLower / Plasma[Conj]Trans
                 */
                else {
                    lda = m < A.lm1 ? A.mb : A.lm%B.mb;
                    CORE_ztrmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(m, m), lda,  /* lda * tempkm */
                               lB,      ldb); /* ldb * tempnn */

                    for (k = m+1; k < A.mt; k++) {
                        tempkm = k == A.mt-1 ? A.m - k * A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);

                        CORE_zgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm,
                            alpha, A(k, m), ldak,
                                   B(k, n), ldbk,
                            zone,  lB,      ldb);
                    }
                }
            }
        }
        /*
         *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
         */
        else {
            if (uplo == PlasmaUpper) {
                if ( trans == PlasmaNoTrans ) {
                    lda = n < A.lm1 ? A.mb : A.lm%B.mb;

                    CORE_ztrmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(n, n), lda,  /* lda * tempkm */
                               lB,      ldb); /* ldb * tempnn */

                    for (k = 0; k < n; k++) {
                        CORE_zgemm(
                            PlasmaNoTrans, trans,
                            tempmm, tempnn, B.mb,
                            alpha, B(m, k), ldb,
                                   A(k, n), A.mb,
                            zone,  lB,      ldb);
                    }
                }
                /*
                 *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
                 */
                else {
                    lda = n < A.lm1 ? A.mb : A.lm%B.mb;

                    CORE_ztrmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(n, n), lda,  /* lda * tempkm */
                               lB,      ldb); /* ldb * tempnn */

                    for (k = n+1; k < A.mt; k++) {
                        tempkn = k == A.nt-1 ? A.n - k * A.nb : A.nb;

                        CORE_zgemm(
                            PlasmaNoTrans, trans,
                            tempmm, tempnn, tempkn,
                            alpha, B(m, k), ldb,
                                   A(n, k), lda,
                            zone,  lB,      ldb);
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaLower / PlasmaNoTrans
             */
            else {
                if (trans == PlasmaNoTrans) {
                    lda = n < A.lm1 ? A.mb : A.lm%B.mb;

                    CORE_ztrmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(n, n), lda,  /* lda * tempkm */
                               lB,      ldb); /* ldb * tempnn */

                    for (k = n+1; k < A.mt; k++) {
                        ldak = BLKLDD(A, k);
                        tempkn = k == A.nt-1 ? A.n - k * A.nb : A.nb;

                        CORE_zgemm(
                            PlasmaNoTrans, trans,
                            tempmm, tempnn, tempkn,
                            alpha, B(m, k), ldb,
                                   A(k, n), ldak,
                            zone,  lB,      ldb);
                    }
                }
                /*
                 *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
                 */
                else {
                    lda = n < A.lm1 ? A.mb : A.lm%B.mb;

                    CORE_ztrmm(
                        side, uplo, trans, diag,
                        tempmm, tempnn,
                        alpha, A(n, n), lda,  /* lda * tempkm */
                               lB,      ldb); /* ldb * tempnn */

                    for (k = 0; k < n; k++) {
                        CORE_zgemm(
                            PlasmaNoTrans, trans,
                            tempmm, tempnn, B.mb,
                            alpha, B(m, k), ldb,
                                   A(n, k), lda,
                            zone,  lB,      ldb);
                    }
                }
            }
        }

        m = next_m;
        n = next_n;
    }
}
#endif

/***************************************************************************//**
 *  Parallel tile triangular matrix-matrix multiplication - dynamic scheduling
 **/
void plasma_pztrmm_quark(PLASMA_enum side, PLASMA_enum uplo,
                         PLASMA_enum trans, PLASMA_enum diag,
                         PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B,
                         PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int lda, ldak, ldb, ldbk;
    int tempkm, tempkn, tempmm, tempnn;

    PLASMA_Complex64_t zone = (PLASMA_Complex64_t)1.0;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    /*
     *  PlasmaLeft / PlasmaUpper / PlasmaNoTrans
     */
    if (side == PlasmaLeft) {
        if (uplo == PlasmaUpper) {
            if (trans == PlasmaNoTrans) {
                for (m = 0; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldb = BLKLDD(B, m);
                    lda = BLKLDD(A, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztrmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(m, m), lda,  /* lda * tempkm */
                                   B(m, n), ldb); /* ldb * tempnn */

                        for (k = m+1; k < A.mt; k++) {
                            tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                            ldbk = BLKLDD(B, k);
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, tempkn, A.mb,
                                alpha, A(m, k), lda,
                                       B(k, n), ldbk,
                                zone,  B(m, n), ldb);
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaUpper / Plasma[Conj]Trans
             */
            else {
                for (m = B.mt-1; m > -1; m--) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldb = BLKLDD(B, m);
                    lda = BLKLDD(A, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztrmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(m, m), lda,  /* lda * tempkm */
                                   B(m, n), ldb); /* ldb * tempnn */

                        for (k = 0; k < m; k++) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, A.mb,
                                alpha, A(k, m), A.mb,
                                       B(k, n), B.mb,
                                zone,  B(m, n), ldb);
                        }
                    }
                }
            }
        }
        /*
         *  PlasmaLeft / PlasmaLower / PlasmaNoTrans
         */
        else {
            if (trans == PlasmaNoTrans) {
                for (m = B.mt-1; m > -1; m--) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldb = BLKLDD(B, m);
                    lda = BLKLDD(A, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztrmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(m, m), lda,  /* lda * tempkm */
                                   B(m, n), ldb); /* ldb * tempnn */

                        for (k = 0; k < m; k++) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, A.mb,
                                alpha, A(m, k), lda,
                                       B(k, n), B.mb,
                                zone,  B(m, n), ldb);
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaLower / Plasma[Conj]Trans
             */
            else {
                for (m = 0; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldb = BLKLDD(B, m);
                    lda = BLKLDD(A, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztrmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(m, m), lda,  /* lda * tempkm */
                                   B(m, n), ldb); /* ldb * tempnn */

                        for (k = m+1; k < A.mt; k++) {
                            tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                            ldak = BLKLDD(A, k);
                            ldbk = BLKLDD(B, k);
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, tempkm, A.mb,
                                alpha, A(k, m), ldak,
                                       B(k, n), ldbk,
                                zone,  B(m, n), ldb);
                        }
                    }
                }
            }
        }
    }
    /*
     *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
     */
    else {
        if (uplo == PlasmaUpper) {
            if (trans == PlasmaNoTrans) {
                for (n = B.nt-1; n > -1; n--) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    lda = BLKLDD(A, n);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldb = BLKLDD(B, m);
                        QUARK_CORE_ztrmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(n, n), lda,  /* lda * tempkm */
                                   B(m, n), ldb); /* ldb * tempnn */

                        for (k = 0; k < n; k++) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, B.mb, A.mb,
                                alpha, B(m, k), ldb,
                                       A(k, n), A.mb,
                                zone,  B(m, n), ldb);
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
             */
            else {
                for (n = 0; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    lda = BLKLDD(A, n);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldb = BLKLDD(B, m);
                        QUARK_CORE_ztrmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(n, n), lda,  /* lda * tempkm */
                                   B(m, n), ldb); /* ldb * tempnn */

                        for (k = n+1; k < A.mt; k++) {
                            tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, tempkn, A.mb,
                                alpha, B(m, k), ldb,
                                       A(n, k), lda,
                                zone,  B(m, n), ldb);
                        }
                    }
                }
            }
        }
        /*
         *  PlasmaRight / PlasmaLower / PlasmaNoTrans
         */
        else {
            if (trans == PlasmaNoTrans) {
                for (n = 0; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    lda = BLKLDD(A, n);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldb = BLKLDD(B, m);
                        QUARK_CORE_ztrmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(n, n), lda,  /* lda * tempkm */
                                   B(m, n), ldb); /* ldb * tempnn */

                        for (k = n+1; k < A.mt; k++) {
                            tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                            ldak = BLKLDD(A, k);
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, tempkn, A.mb,
                                alpha, B(m, k), ldb,
                                       A(k, n), ldak,
                                zone,  B(m, n), ldb);
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
             */
            else {
                for (n = B.nt-1; n > -1; n--) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    lda = BLKLDD(A, n);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldb = BLKLDD(B, m);
                        QUARK_CORE_ztrmm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A.mb,
                            alpha, A(n, n), lda,  /* lda * tempkm */
                                   B(m, n), ldb); /* ldb * tempnn */

                        for (k = 0; k < n; k++) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, B.mb, A.mb,
                                alpha, B(m, k), ldb,
                                       A(n, k), lda,
                                zone,  B(m, n), ldb);
                        }
                    }
                }
            }
        }
    }
}
