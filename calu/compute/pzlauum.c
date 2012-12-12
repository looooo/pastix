/**
 *
 * @file pzlauum.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel UU' or L'L operation - dynamic scheduling
 **/
void plasma_pzlauum_quark(PLASMA_enum uplo, PLASMA_desc A,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldam;
    int tempkm, tempmm, tempnn;

    PLASMA_Complex64_t zone = (PLASMA_Complex64_t)1.0;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    /*
     *  PlasmaLower
     */
    if (uplo == PlasmaLower) {
        for (m = 0; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            for(n = 0; n < m; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_zherk(
                    plasma->quark, &task_flags,
                    uplo, PlasmaConjTrans,
                    tempnn, tempmm, A.mb,
                    1.0, A(m, n), ldam,
                    1.0, A(n, n), A.mb);

                for(k = n+1; k < m; k++) {
                    tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                    QUARK_CORE_zgemm(
                        plasma->quark, &task_flags,
                        PlasmaConjTrans, PlasmaNoTrans,
                        tempkm, tempnn, tempmm, A.mb,
                        zone, A(m, k), ldam,
                              A(m, n), ldam,
                        zone, A(k, n), A.mb);
                }
            }
            for (n = 0; n < m; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_ztrmm(
                    plasma->quark, &task_flags,
                    PlasmaLeft, uplo, PlasmaConjTrans, PlasmaNonUnit,
                    tempmm, tempnn, A.mb,
                    zone, A(m, m), ldam,
                          A(m, n), ldam);
            }
            QUARK_CORE_zlauum(
                plasma->quark, &task_flags,
                uplo,
                tempmm,
                A.mb, A(m, m), ldam);
        }
    }
    /*
     *  PlasmaUpper
     */
    else {
        for (m = 0; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            for (n = 0; n < m; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_zherk(
                    plasma->quark, &task_flags,
                    uplo, PlasmaNoTrans,
                    tempnn, tempmm, A.mb,
                    1.0, A(n, m), A.mb,
                    1.0, A(n, n), A.mb);

                for (k = n+1; k < m; k++){
                    tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                    QUARK_CORE_zgemm(
                        plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaConjTrans,
                        tempnn, tempkm, tempmm, A.mb,
                        zone, A(n, m), A.mb,
                              A(k, m), A.mb,
                        zone, A(n, k), A.mb);
                }
            }
            for (n = 0; n < m; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_ztrmm(
                    plasma->quark, &task_flags,
                    PlasmaRight, uplo, PlasmaConjTrans, PlasmaNonUnit,
                    tempnn, tempmm, A.mb,
                    zone, A(m, m), ldam,
                          A(n, m), A.mb);
            }
            QUARK_CORE_zlauum(
                plasma->quark, &task_flags,
                uplo,
                tempmm,
                A.mb, A(m, m), ldam);
        }
    }
}
