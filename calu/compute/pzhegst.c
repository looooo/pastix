/**
 *
 * @file pzhegst.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define B(m,n) BLKADDR(B, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel Transformation to standard eigenvalue problem  - dynamic scheduler
 **/
void plasma_pzhegst_quark(PLASMA_enum itype, PLASMA_enum uplo,
                          PLASMA_desc A, PLASMA_desc B,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k;
    int ldak, ldbk;
    int tempkn;
    static double done = 1.0;
    static PLASMA_Complex64_t zone = 1.0;
    static PLASMA_Complex64_t mzone = -1.0;
    static PLASMA_Complex64_t zhalf = 0.5;
    static PLASMA_Complex64_t mzhalf = -0.5;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;

    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    if (itype == 1) {
        if (uplo == PlasmaLower) {
            /*
             * Type 1 / PlasmaLower
             */
            for (k = 0; k < A.nt; k++){
                tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);

                QUARK_CORE_zhegst(
                    plasma->quark, &task_flags,
                    itype, uplo, tempkn,
                    A(k, k), ldak,
                    B(k, k), ldbk,
                    sequence, request, A.nb*k);

                if (k*A.nb+tempkn < A.n) {
                    plasma_pztrsm_quark(
                        PlasmaRight, uplo,
                        PlasmaConjTrans, PlasmaNonUnit,
                        zone,
                        plasma_desc_submatrix(B, k*B.nb,        k*B.nb, tempkn,            tempkn),
                        plasma_desc_submatrix(A, k*A.nb+tempkn, k*A.nb, A.n-k*A.nb-tempkn, tempkn),
                        sequence, request);

                    plasma_pzhemm_quark(
                        PlasmaRight, uplo, mzhalf,
                        plasma_desc_submatrix(A, k*A.nb,        k*A.nb, tempkn,            tempkn),
                        plasma_desc_submatrix(B, k*B.nb+tempkn, k*B.nb, B.n-k*B.nb-tempkn, tempkn),
                        zone,
                        plasma_desc_submatrix(A, k*A.nb+tempkn, k*A.nb, A.n-k*A.nb-tempkn, tempkn),
                        sequence, request);

                    plasma_pzher2k_quark(
                        uplo, PlasmaNoTrans,
                        mzone,
                        plasma_desc_submatrix(A, k*A.nb+tempkn, k*A.nb,        A.n-k*A.nb-tempkn, tempkn),
                        plasma_desc_submatrix(B, k*B.nb+tempkn, k*B.nb,        B.n-k*B.nb-tempkn, tempkn),
                        done,
                        plasma_desc_submatrix(A, k*A.nb+tempkn, k*A.nb+tempkn, A.n-k*A.nb-tempkn, A.n-k*A.nb-tempkn),
                        sequence, request);

                    plasma_pzhemm_quark(
                        PlasmaRight, uplo,
                        mzhalf,
                        plasma_desc_submatrix(A, k*A.nb,        k*A.nb, tempkn,            tempkn),
                        plasma_desc_submatrix(B, k*B.nb+tempkn, k*B.nb, B.n-k*B.nb-tempkn, tempkn),
                        zone,
                        plasma_desc_submatrix(A, k*A.nb+tempkn, k*A.nb, A.n-k*A.nb-tempkn, tempkn),
                        sequence, request);

                    plasma_pztrsm_quark(
                        PlasmaLeft, uplo, PlasmaNoTrans, PlasmaNonUnit, zone,
                        plasma_desc_submatrix(B, k*B.nb+tempkn, k*B.nb+tempkn, B.n-k*B.nb-tempkn, B.n-k*B.nb-tempkn),
                        plasma_desc_submatrix(A, k*A.nb+tempkn, k*A.nb,        A.n-k*A.nb-tempkn, tempkn),
                        sequence, request);
                }
            }
        }
        else {
            /*
             * Type 1 / PlasmaUpper
             */
            for (k = 0; k < A.nt; k++){
                tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);
                QUARK_CORE_zhegst(
                    plasma->quark, &task_flags,
                    itype, uplo, tempkn,
                    A(k, k), ldak,
                    B(k, k), ldbk,
                    sequence, request, A.nb*k);

                if (k*A.nb+tempkn < A.n) {
                    plasma_pztrsm_quark(
                        PlasmaLeft, uplo,
                        PlasmaConjTrans, PlasmaNonUnit,
                        zone,
                        plasma_desc_submatrix(B, k*B.nb, k*B.nb,        tempkn, tempkn),
                        plasma_desc_submatrix(A, k*A.nb, k*A.nb+tempkn, tempkn, A.n-k*A.nb-tempkn),
                        sequence, request);

                    plasma_pzhemm_quark(
                        PlasmaLeft, uplo, mzhalf,
                        plasma_desc_submatrix(A, k*A.nb, k*A.nb,        tempkn, tempkn),
                        plasma_desc_submatrix(B, k*B.nb, k*B.nb+tempkn, tempkn, B.n-k*B.nb-tempkn),
                        zone,
                        plasma_desc_submatrix(A, k*A.nb, k*A.nb+tempkn, tempkn, A.n-k*A.nb-tempkn),
                        sequence, request);

                    plasma_pzher2k_quark(
                        uplo, PlasmaConjTrans,
                        mzone,
                        plasma_desc_submatrix(A, k*A.nb,        k*A.nb+tempkn, tempkn,            A.n-k*A.nb-tempkn),
                        plasma_desc_submatrix(B, k*B.nb,        k*B.nb+tempkn, tempkn,            B.n-k*B.nb-tempkn),
                        done,
                        plasma_desc_submatrix(A, k*A.nb+tempkn, k*A.nb+tempkn, A.n-k*A.nb-tempkn, A.n-k*A.nb-tempkn),
                        sequence, request);

                    plasma_pzhemm_quark(
                        PlasmaLeft, uplo,
                        mzhalf,
                        plasma_desc_submatrix(A, k*A.nb, k*A.nb,        tempkn, tempkn),
                        plasma_desc_submatrix(B, k*B.nb, k*B.nb+tempkn, tempkn, B.n-k*B.nb-tempkn),
                        zone,
                        plasma_desc_submatrix(A, k*A.nb, k*A.nb+tempkn, tempkn, A.n-k*A.nb-tempkn),
                        sequence, request);

                    plasma_pztrsm_quark(
                        PlasmaRight, uplo,
                        PlasmaNoTrans, PlasmaNonUnit,
                        zone,
                        plasma_desc_submatrix(B, k*B.nb+tempkn, k*B.nb+tempkn, A.n-k*A.nb-tempkn, A.n-k*A.nb-tempkn),
                        plasma_desc_submatrix(A, k*A.nb,        k*A.nb+tempkn, tempkn,            A.n-k*A.nb-tempkn),
                        sequence, request);
                }
            }
        }
    }
    else{
        if (uplo == PlasmaLower) {
            /*
             * Type 2 or 3  / PlasmaLower
             */
            for (k = 0; k < A.nt; k++){
                tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);

                plasma_pztrmm_quark(
                    PlasmaRight, uplo,
                    PlasmaNoTrans, PlasmaNonUnit,
                    zone,
                    plasma_desc_submatrix(B, 0,      0, k*B.nb, k*B.nb),
                    plasma_desc_submatrix(A, k*A.nb, 0, tempkn, k*A.nb),
                    sequence, request);

                plasma_pzhemm_quark(
                    PlasmaLeft, uplo, zhalf,
                    plasma_desc_submatrix(A, k*A.nb, k*A.nb, tempkn, tempkn),
                    plasma_desc_submatrix(B, k*B.nb, 0,      tempkn, k*B.nb),
                    zone,
                    plasma_desc_submatrix(A, k*A.nb, 0,      tempkn, k*A.nb),
                    sequence, request);

                plasma_pzher2k_quark(
                    uplo, PlasmaConjTrans,
                    zone,
                    plasma_desc_submatrix(A, k*A.nb, 0, tempkn, k*A.nb),
                    plasma_desc_submatrix(B, k*B.nb, 0, tempkn, k*B.nb),
                    done,
                    plasma_desc_submatrix(A, 0, 0, k*A.nb, k*A.nb),
                    sequence, request);

                plasma_pzhemm_quark(
                    PlasmaLeft, uplo, zhalf,
                    plasma_desc_submatrix(A, k*A.nb, k*A.nb, tempkn, tempkn),
                    plasma_desc_submatrix(B, k*B.nb, 0,      tempkn, k*B.nb),
                    zone,
                    plasma_desc_submatrix(A, k*A.nb, 0,      tempkn, k*A.nb),
                    sequence, request);

                plasma_pztrmm_quark(
                    PlasmaLeft, uplo,
                    PlasmaConjTrans, PlasmaNonUnit,
                    zone,
                    plasma_desc_submatrix(B, k*B.nb, k*B.nb, tempkn, tempkn),
                    plasma_desc_submatrix(A, k*A.nb, 0,      tempkn, k*A.nb),
                    sequence, request);

                QUARK_CORE_zhegst(
                    plasma->quark, &task_flags,
                    itype, uplo, tempkn,
                    A(k, k), ldak,
                    B(k, k), ldbk,
                    sequence, request, A.nb*k);
            }
        }
        else {
            /*
             * Type 2 or 3  / PlasmaUpper
             */
            for (k = 0; k < A.nt; k++){
                tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);

                plasma_pztrmm_quark(
                    PlasmaLeft, uplo,
                    PlasmaNoTrans, PlasmaNonUnit,
                    zone,
                    plasma_desc_submatrix(B, 0, 0,      k*B.nb, k*B.nb),
                    plasma_desc_submatrix(A, 0, k*A.nb, k*A.nb, tempkn),
                    sequence, request);

                plasma_pzhemm_quark(
                    PlasmaRight, uplo, zhalf,
                    plasma_desc_submatrix(A, k*A.nb, k*A.nb, tempkn, tempkn),
                    plasma_desc_submatrix(B, 0,      k*B.nb, k*B.nb, tempkn),
                    zone,
                    plasma_desc_submatrix(A, 0,      k*A.nb, k*A.nb, tempkn),
                    sequence, request);

                plasma_pzher2k_quark(
                    uplo, PlasmaNoTrans, zone,
                    plasma_desc_submatrix(A, 0, k*A.nb, k*A.nb, tempkn),
                    plasma_desc_submatrix(B, 0, k*B.nb, k*B.nb, tempkn),
                    done,
                    plasma_desc_submatrix(A, 0, 0,      k*A.nb, k*A.nb),
                    sequence, request);

                plasma_pzhemm_quark(
                    PlasmaRight, uplo, zhalf,
                    plasma_desc_submatrix(A, k*A.nb, k*A.nb, tempkn, tempkn),
                    plasma_desc_submatrix(B, 0,      k*B.nb, k*B.nb, tempkn),
                    zone,
                    plasma_desc_submatrix(A, 0,      k*A.nb, k*A.nb, tempkn),
                    sequence, request);

                plasma_pztrmm_quark(
                    PlasmaRight, uplo,
                    PlasmaConjTrans, PlasmaNonUnit,
                    zone,
                    plasma_desc_submatrix(B, k*B.nb, k*B.nb, tempkn, tempkn),
                    plasma_desc_submatrix(A, 0,      k*A.nb, k*A.nb, tempkn),
                    sequence, request);

                QUARK_CORE_zhegst(
                    plasma->quark, &task_flags,
                    itype, uplo, tempkn,
                    A(k, k), ldak,
                    B(k, k), ldbk,
                    sequence, request, A.nb*k);
            }
        }
    }
}
