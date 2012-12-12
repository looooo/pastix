/**
 *
 * @file pztrsm.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define B(m,n) BLKADDR(B, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel tile triangular solve - static scheduling
 **/
void plasma_pztrsm(plasma_context_t *plasma)
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
    int next_k;
    int next_m;
    int next_n;
    int lda, ldb;
    int tempkm, tempnn, tempmm, tempkn;

    PLASMA_Complex64_t zone  = (PLASMA_Complex64_t) 1.0;
    PLASMA_Complex64_t mzone = (PLASMA_Complex64_t)-1.0;
    PLASMA_Complex64_t lalpha;
    PLASMA_Complex64_t minvalpha;

    plasma_unpack_args_9(side, uplo, trans, diag, alpha, A, B, sequence, request);
    minvalpha = mzone / alpha;
    if (sequence->status != PLASMA_SUCCESS)
        return;
    ss_init(B.mt, B.nt, -1);
    /*
     *  PlasmaLeft
     */
    if (side == PlasmaLeft) {
        k = 0;
        m = PLASMA_RANK;
        while (m >= B.mt) {
            k++;
            m = m - B.mt + k;
        }
        n = 0;

        while (k < B.mt && m < B.mt) {
            next_n = n;
            next_m = m;
            next_k = k;

            next_n++;
            if (next_n >= B.nt) {
                next_m += PLASMA_SIZE;
                while (next_m >= B.mt && next_k < B.mt) {
                    next_k++;
                    next_m = next_m - B.mt + next_k;
                }
                next_n = 0;
            }

            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
            tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;

            lalpha = k == 0 ? alpha : zone;
            if (m == k) {
                ss_cond_wait(m, n, k-1);
                /*
                 *  PlasmaLeft / PlasmaLower / PlasmaNoTrans
                 *  PlasmaLeft / PlasmaUpper / Plasma[Conj]Trans
                 */
                if ((uplo == PlasmaLower && trans == PlasmaNoTrans)
                    || (uplo == PlasmaUpper && trans != PlasmaNoTrans)) {
                    tempkm = k == B.mt-1 ? B.m-k*B.mb : B.mb;
                    lda = BLKLDD(A, k);
                    ldb = BLKLDD(B, k);
                    CORE_ztrsm(
                        side, uplo, trans, diag,
                        tempkm, tempnn,
                        lalpha, A(k, k), lda,
                                B(k, n), ldb);
                }
                /*
                 *  PlasmaLeft / PlasmaLower / Plasma[Cojn]Trans
                 *  PlasmaLeft / PlasmaUpper / PlasmaNoTrans
                 */
                else {
                    tempkm = k == 0 ? B.m-(B.mt-1)*B.mb : B.mb;
                    lda = BLKLDD(A, B.mt-1-k);
                    ldb = BLKLDD(B, B.mt-1-k);
                    CORE_ztrsm(
                        side, uplo, trans, diag,
                        tempkm, tempnn,
                        lalpha, A(B.mt-1-k, B.mt-1-k), lda,
                                B(B.mt-1-k, n       ), ldb);
                }
                ss_cond_set(k, n, k);
            }
            else {
                ss_cond_wait(k, n, k);
                ss_cond_wait(m, n, k-1);
                /*
                 *  PlasmaRight / PlasmaLower / PlasmaNoTrans
                 */
                if (uplo == PlasmaLower) {
                    if (trans == PlasmaNoTrans) {
                        lda = BLKLDD(A, m);
                        ldb = BLKLDD(B, m);
                        CORE_zgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            tempmm, tempnn, B.mb,
                            mzone,  A(m, k), lda,
                                    B(k, n), B.mb,
                            lalpha, B(m, n), ldb);
                    }
                    /*
                     *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
                     */
                    else {
                        tempkm = k == 0 ? A.m-(A.mt-1)*A.mb : A.mb;
                        lda = BLKLDD(A, B.mt-1-k);
                        ldb = BLKLDD(B, B.mt-1-k);
                        CORE_zgemm(
                            trans, PlasmaNoTrans,
                            B.mb, tempnn, tempkm,
                            mzone,  A(A.mt-1-k, A.mt-1-m), lda,
                                    B(B.mt-1-k, n       ), ldb,
                            lalpha, B(B.mt-1-m, n       ), B.mb);
                    }
                }
                else {
                    /*
                     *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
                     */
                    if (trans == PlasmaNoTrans) {
                        tempkm = k == 0 ? A.m-(A.mt-1)*A.mb : A.mb;
                        ldb = BLKLDD(B, B.mt-1-k);
                        CORE_zgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            B.mb, tempnn, tempkm,
                            mzone,  A(A.mt-1-m, A.mt-1-k), A.mb,
                                    B(B.mt-1-k, n       ), ldb,
                            lalpha, B(B.mt-1-m, n       ), B.mb);
                    }
                    /*
                     *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
                     */
                    else {
                        ldb = BLKLDD(B, m);
                        CORE_zgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, B.mb,
                            mzone,  A(k, m), A.mb,
                                    B(k, n), B.mb,
                            lalpha, B(m, n), ldb);
                    }
                }
                ss_cond_set(m, n, k);
            }
            n = next_n;
            m = next_m;
            k = next_k;
        }
    }
    /*
     *  PlasmaRight
     */
    else {
        k = 0;
        n = PLASMA_RANK;
        while (n >= B.nt) {
            k++;
            n = n - B.nt + k;
        }
        m = 0;

        while (k < B.nt && n < B.nt) {
            next_n = n;
            next_m = m;
            next_k = k;

            next_m++;
            if (next_m >= B.mt) {
                next_n += PLASMA_SIZE;
                while (next_n >= B.nt && next_k < B.nt) {
                    next_k++;
                    next_n = next_n - B.nt + next_k;
                }
                next_m = 0;
            }

            tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;

            lalpha = k == 0 ? alpha : zone;
            if (n == k) {
                ss_cond_wait(m, n, k-1);
                /*
                 *  PlasmaRight / PlasmaLower / PlasmaNoTrans
                 */
                if (uplo == PlasmaLower) {
                    if (trans == PlasmaNoTrans) {
                        tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                        lda = BLKLDD(A, B.nt-1-k);
                        ldb = BLKLDD(B, m);
                        CORE_ztrsm(
                            side, uplo, trans, diag,
                            tempmm, tempkn,
                            lalpha, A(B.nt-1-k, B.nt-1-k), lda,
                                    B(m,        B.nt-1-k), ldb);
                    }
                    /*
                     *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
                     */
                    else {
                        tempkn = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                        lda = BLKLDD(A, k);
                        ldb = BLKLDD(B, m);
                        CORE_ztrsm(
                            side, uplo, trans, diag,
                            tempmm, tempkn,
                            alpha, A(k, k), lda,
                                   B(m, k), ldb);
                    }
                }
                else {
                    /*
                     *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
                     */
                    if (trans == PlasmaNoTrans) {
                        tempkn = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                        lda = BLKLDD(A, k);
                        ldb = BLKLDD(B, m);
                        CORE_ztrsm(
                            side, uplo, trans, diag,
                            tempmm, tempkn,
                            lalpha, A(k, k), lda,
                                    B(m, k), ldb);
                    }
                    /*
                     *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
                     */
                    else {
                        tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                        lda = BLKLDD(A, B.nt-1-k);
                        ldb = BLKLDD(B, m);
                        CORE_ztrsm(
                            side, uplo, trans, diag,
                            tempmm, tempkn,
                            alpha, A(B.nt-1-k, B.nt-1-k), lda,
                                   B(m,        B.nt-1-k), ldb);
                    }
                }
                ss_cond_set(m, k, k);
            }
            else {
                ss_cond_wait(m, k, k);
                ss_cond_wait(m, n, k-1);
                /*
                 *  PlasmaRight / PlasmaLower / PlasmaNoTrans
                 */
                if (uplo == PlasmaLower) {
                    if (trans == PlasmaNoTrans) {
                        tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                        lda = BLKLDD(A, B.nt-1-k);
                        ldb = BLKLDD(B, m);
                        CORE_zgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            tempmm, B.mb, tempkn,
                            mzone,  B(m,        B.nt-1-k), ldb,
                                    A(B.nt-1-k, B.nt-1-n), lda,
                            lalpha, B(m,        B.nt-1-n), ldb);
                    }
                    /*
                     *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
                     */
                    else {
                        lda = BLKLDD(A, n);
                        ldb = BLKLDD(B, m);
                        CORE_zgemm(
                            PlasmaNoTrans, trans,
                            tempmm, tempnn, B.mb,
                            minvalpha, B(m, k), ldb,
                                       A(n, k), lda,
                            zone,      B(m, n), ldb);
                    }
                }
                else {
                    /*
                     *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
                     */
                    if (trans == PlasmaNoTrans) {
                        lda = BLKLDD(A, k);
                        ldb = BLKLDD(B, m);
                        CORE_zgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            tempmm, tempnn, B.mb,
                            mzone,  B(m, k), ldb,
                                    A(k, n), lda,
                            lalpha, B(m, n), ldb);
                    }
                    /*
                     *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
                     */
                    else {
                        tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                        ldb = BLKLDD(B, m);
                        CORE_zgemm(
                            PlasmaNoTrans, trans,
                            tempmm, B.nb, tempkn,
                            minvalpha, B(m,        B.nt-1-k), ldb,
                                       A(B.nt-1-n, B.nt-1-k), A.mb,
                            zone,      B(m,        B.nt-1-n), ldb);
                    }
                }
                ss_cond_set(m, n, k);
            }
            n = next_n;
            m = next_m;
            k = next_k;
        }
    }
    ss_finalize();
}

/***************************************************************************//**
 *  Parallel tile triangular solve - dynamic scheduling
 **/
void plasma_pztrsm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum trans, PLASMA_enum diag,
                         PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B,
                         PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int lda, ldan, ldb;
    int tempkm, tempkn, tempmm, tempnn;

    PLASMA_Complex64_t zone       = (PLASMA_Complex64_t) 1.0;
    PLASMA_Complex64_t mzone      = (PLASMA_Complex64_t)-1.0;
    PLASMA_Complex64_t minvalpha  = (PLASMA_Complex64_t)-1.0 / alpha;
    PLASMA_Complex64_t lalpha;

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
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == 0 ? B.m-(B.mt-1)*B.mb : B.mb;
                    lda = BLKLDD(A, B.mt-1-k);
                    ldb = BLKLDD(B, B.mt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztrsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A.mb,
                            lalpha, A(B.mt-1-k, B.mt-1-k), lda,  /* lda * tempkm */
                                    B(B.mt-1-k,        n), ldb); /* ldb * tempnn */
                    }
                    for (m = k+1; m < B.mt; m++) {
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaNoTrans,
                                B.mb, tempnn, tempkm, A.mb,
                                mzone,  A(B.mt-1-m, B.mt-1-k), A.mb,
                                        B(B.mt-1-k, n       ), ldb,
                                lalpha, B(B.mt-1-m, n       ), B.mb);
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaUpper / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == B.mt-1 ? B.m-k*B.mb : B.mb;
                    lda = BLKLDD(A, k);
                    ldb = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztrsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A.mb,
                            lalpha, A(k, k), lda,
                                    B(k, n), ldb);
                    }
                    for (m = k+1; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldb = BLKLDD(B, m);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, A.mb,
                                mzone,  A(k, m), A.mb,
                                        B(k, n), B.mb,
                                lalpha, B(m, n), ldb);
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
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == B.mt-1 ? B.m-k*B.mb : B.mb;
                    lda = BLKLDD(A, k);
                    ldb = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztrsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A.mb,
                            lalpha, A(k, k), lda,
                                    B(k, n), ldb);
                    }
                    for (m = k+1; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        lda = BLKLDD(A, m);
                        ldb = BLKLDD(B, m);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, A.mb,
                                mzone,  A(m, k), lda,
                                        B(k, n), B.mb,
                                lalpha, B(m, n), ldb);
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaLower / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == 0 ? B.m-(B.mt-1)*B.mb : B.mb;
                    lda = BLKLDD(A, B.mt-1-k);
                    ldb = BLKLDD(B, B.mt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztrsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A.mb,
                            lalpha, A(B.mt-1-k, B.mt-1-k), lda,
                                    B(B.mt-1-k,        n), ldb);
                    }
                    for (m = k+1; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                B.mb, tempnn, tempkm, A.mb,
                                mzone,  A(B.mt-1-k, B.mt-1-m), lda,
                                        B(B.mt-1-k, n       ), ldb,
                                lalpha, B(B.mt-1-m, n       ), B.mb);
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
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                    lda = BLKLDD(A, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldb = BLKLDD(B, m);
                        QUARK_CORE_ztrsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A.mb,
                            lalpha, A(k, k), lda,  /* lda * tempkn */
                                    B(m, k), ldb); /* ldb * tempkn */
                    }
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldb = BLKLDD(B, m);
                        for (n = k+1; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, A.mb,
                                mzone,  B(m, k), ldb,  /* ldb * B.mb   */
                                        A(k, n), lda,  /* lda * tempnn */
                                lalpha, B(m, n), ldb); /* ldb * tempnn */
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                    lda = BLKLDD(A, B.nt-1-k);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldb = BLKLDD(B, m);
                        QUARK_CORE_ztrsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A.mb,
                            alpha, A(B.nt-1-k, B.nt-1-k), lda,  /* lda * tempkn */
                                   B(       m, B.nt-1-k), ldb); /* ldb * tempkn */

                        for (n = k+1; n < B.nt; n++) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, B.nb, tempkn, A.mb,
                                minvalpha, B(m,        B.nt-1-k), ldb,  /* ldb  * tempkn */
                                           A(B.nt-1-n, B.nt-1-k), A.mb, /* A.mb * tempkn (Never last row) */
                                zone,      B(m,        B.nt-1-n), ldb); /* ldb  * B.nb   */
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
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                    lda = BLKLDD(A, B.nt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldb = BLKLDD(B, m);
                        QUARK_CORE_ztrsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A.mb,
                            lalpha, A(B.nt-1-k, B.nt-1-k), lda,  /* lda * tempkn */
                                    B(       m, B.nt-1-k), ldb); /* ldb * tempkn */

                        for (n = k+1; n < B.nt; n++) {
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, B.nb, tempkn, A.mb,
                                mzone,  B(m,        B.nt-1-k), ldb,  /* ldb * tempkn */
                                        A(B.nt-1-k, B.nt-1-n), lda,  /* lda * B.nb   */
                                lalpha, B(m,        B.nt-1-n), ldb); /* ldb * B.nb   */
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                    lda = BLKLDD(A, k);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldb = BLKLDD(B, m);
                        QUARK_CORE_ztrsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A.mb,
                            alpha, A(k, k), lda,  /* lda * tempkn */
                                   B(m, k), ldb); /* ldb * tempkn */

                        for (n = k+1; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            ldan = BLKLDD(A, n);
                            QUARK_CORE_zgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, B.mb, A.mb,
                                minvalpha, B(m, k), ldb,  /* ldb  * tempkn */
                                           A(n, k), ldan, /* ldan * tempkn */
                                zone,      B(m, n), ldb); /* ldb  * tempnn */
                        }
                    }
                }
            }
        }
    }
}
