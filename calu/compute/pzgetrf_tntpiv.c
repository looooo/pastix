/**
 *
 * @file pzgetrf_reclap.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * LU with Partial pivoting.
 *
 * @version 2.4.5
 * @author Mathieu Faverge
 * @author Hatem Ltaief
 * @date 2009-11-15
 *
 * @precisions normal z -> s d c
 *
 **/
#include <lapacke.h>
#include "common.h"

#define A(m,n)    (BLKADDR(A, PLASMA_Complex64_t, m, n))
#define Acpy(__m) (BLKADDR(W, PLASMA_Complex64_t, (__m), 0))
#define IPIV(k)    (&(IPIV[(int64_t)A.mb*(int64_t)(k)]))
#define RANK(__m, __k) (&(Wi[(int64_t)W.mb*(int64_t)(__m) + (int64_t)W.lm*(__k)]))

void plasma_pzgetrf_rectil_panel_quark(plasma_context_t *plasma,
                                       int *panel_thread_count,
                                       PLASMA_desc A, int *IPIV,
                                       Quark_Task_Flags *task_flags,
                                       PLASMA_sequence *sequence, PLASMA_request *request)
{
    while ( ((*panel_thread_count * 4 * A.mb) > A.m)
            && (*panel_thread_count > 1) ) {
        *panel_thread_count--;
        QUARK_Task_Flag_Set(task_flags, TASK_THREAD_COUNT, *panel_thread_count );
    }

    QUARK_CORE_zgetrf_rectil(
        plasma->quark, task_flags,
        A, A(0, 0), A.mb*A.nb, IPIV,
        sequence, request, 1, A.i,
        *panel_thread_count );
}

void plasma_pzgetrf_tntpiv_panel_quark(plasma_context_t *plasma,
                                       PLASMA_desc A, int *IPIV,
                                       PLASMA_desc W, int *Wi,
                                       Quark_Task_Flags *task_flags,
                                       PLASMA_sequence *sequence, PLASMA_request *request)
{
    int tempkm, tempmm, tempnn, tempr;
    int tempm, nexti, pos;
    int ldak, ldam;
    int round, round_size = 4;
    int prev_round_size = 1;
    int curr_round_size = 4;
    int next_round_size = curr_round_size * round_size;
    int m, i;

    tempkm = min(A.m, A.mb);
    ldak = BLKLDD(A, 0);

    /* Create a first copy of the panel */
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m - m*A.mb : A.mb;
        ldam = BLKLDD(A, m);
        i   = m / round_size;
        pos = m % round_size;

        QUARK_CORE_zlacpy_f1(
            plasma->quark, task_flags,
            PlasmaUpperLower,
            tempmm, A.n, A.mb,
            A(m, 0), ldam,
            Acpy(i) + pos*A.mb, W.mb,
            Acpy(i), W.mb*W.nb, OUTPUT | GATHERV );
    }

    round = 0;
    while ( ((A.mt - 1) / curr_round_size) > 0 ) {
        /* Let's submit all the factorizations */
        for (m=0, i=0; m<A.mt; m+=curr_round_size, i++) {

            if ( (m+curr_round_size) < A.mt ) {
                tempm = round_size * A.mb;
                tempr = curr_round_size * A.mb;
            } else {
                tempm = ((A.mt-1-m) / prev_round_size) * A.mb;
                if ( ( (A.mt-1-m) % prev_round_size == 0 )
                     && (A.m % A.mb != 0) ) {
                    tempm += A.m % A.mb;
                } else {
                    tempm += A.mb;
                }
                tempr = A.m - m * A.mb;
            }

            QUARK_CORE_zgetrf(
                plasma->quark, task_flags,
                tempm, A.n, A.mb,
                Acpy( i ), W.mb,
                IPIV( i ),
                sequence, request,
                0, A.i );

            tempm = min(tempm, A.n);
            nexti = i / round_size;
            pos   = i % round_size;

            if ( pos == 0 ) {
                assert( nexti <= i );
                QUARK_CORE_zlaset(
                    plasma->quark, task_flags,
                    PlasmaUpperLower, W.mb, W.nb, 0., 0., Acpy( nexti ), W.mb );
            }

            QUARK_CORE_zlacpy_pivot(
                plasma->quark, task_flags,
                plasma_desc_submatrix(A, m*A.mb, 0, tempr, A.n),
                1, tempm, IPIV( i ),
                RANK( i, round ), RANK( nexti, round+1 ),
                Acpy( nexti ), W.mb,
                pos*A.mb, prev_round_size==1 );
        }

        round++;
        prev_round_size = curr_round_size;
        curr_round_size = next_round_size;
        next_round_size = curr_round_size * round_size;
    }

    /* Last factorization */
    tempm = ((A.mt-1) / prev_round_size) * A.mb;
    if ( ( (A.mt-1) % prev_round_size == 0 )
         && (A.m % A.mb != 0) )
        tempm += A.m % A.mb;
    else
        tempm += A.mb;

    QUARK_CORE_zgetrf(
        plasma->quark, task_flags,
        tempm, A.n, A.mb,
        Acpy(0), W.mb,
        IPIV(0),
        sequence, request,
        1, A.i );

    QUARK_CORE_pivot_update(
        plasma->quark, task_flags,
        tempm, min(tempm, A.n),
        RANK( 0, round ), IPIV( 0 ),
        A.i, (int)(prev_round_size == 1));

    /* Finish to factorize the panel */
    QUARK_CORE_zlaswp_ontile(
        plasma->quark, task_flags,
        A, A(0, 0), 1, min(tempkm, A.n), IPIV(0), 1, A(0, 0) );

    /* Copy back the factorization result, once it has been swapped */
    QUARK_CORE_zlacpy(
        plasma->quark, task_flags,
        PlasmaUpperLower,
        tempkm, A.n, A.mb,
        Acpy(0), W.mb,
        A(0, 0), ldak);

    /* Apply TRSM on the panel
     * Using A(k,k) ensures that the panel is swapped */
    for (m=1; m<A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        ldam = BLKLDD(A, m);

        QUARK_CORE_ztrsm(
            plasma->quark, task_flags,
            PlasmaRight, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit,
            tempmm, A.n, A.mb,
            1., A(0, 0), ldak,
                A(m, 0), ldam);
    }
}

/*
 * W is a workspace in tile layout with tile of size max_round_size*A.mb -by- A.nb
 * Wi is an integer workspace to store the rank of the lines involved in each round.
 * This workspace has tiles of size max_round_size*A.mb -by- 1
 */

/***************************************************************************//**
 *  Parallel tile LU factorization - dynamic scheduling - Right looking
 **/
void plasma_pzgetrf_tntpiv_quark(PLASMA_desc A, int *IPIV,
                                 PLASMA_desc W, int *Wi,
                                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    int i, k, m, n, minmnt;
    plasma_context_t *plasma;
    int tempkm, tempkn, tempmm, tempnn;
    int tempm, tempk;
    int ldak, ldam;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    PLASMA_Complex64_t zone  = (PLASMA_Complex64_t)1.0;
    PLASMA_Complex64_t mzone = (PLASMA_Complex64_t)-1.0;
    void * fakedep;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    int ok = 1;
    PLASMA_Complex64_t *A0 = NULL;
    A0 = (PLASMA_Complex64_t*)malloc( (A.m) * (A.n) * sizeof(PLASMA_Complex64_t) );
    if ( ! A0 ) {
        fprintf(stderr, "Our of Memory for %s\n", "A0");
        return;
    }
    PLASMA_zplrnt(A.m, A.n, A0, A.m, 4361);

    minmnt = min(A.mt, A.nt);
    for (k = 0; k < minmnt; k++)
    {
        tempk  = k * A.mb;
        tempm  = A.m - tempk;
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        ldak = BLKLDD(A, k);

        QUARK_Task_Flag_Set(&task_flags, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - k );

        plasma_pzgetrf_tntpiv_panel_quark(plasma,
                                          plasma_desc_submatrix(A, tempk, tempk, tempm, tempkn),
                                          IPIV(k), W, Wi, &task_flags,
                                          sequence, request);

        /*
         * Update the trailing submatrix
         */
        fakedep = (void *)(intptr_t)(k+1);
        for (n = k+1; n < A.nt; n++)
        {

            QUARK_Task_Flag_Set(&task_flags, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - n );
            /*
             * Apply row interchange after the panel (work on the panel)
             */
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;

            QUARK_CORE_zswptr_ontile(
                plasma->quark, &task_flags,
                plasma_desc_submatrix(A, tempk, n*A.nb, tempm, tempnn),
                A(k, n), 1, tempkm, IPIV(k), 1,
                A(k, k), ldak);

            m = k+1;
            if ( m < A.mt ) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);

                QUARK_CORE_zgemm2(
                    plasma->quark, &task_flags,
                    PlasmaNoTrans, PlasmaNoTrans,
                    tempmm, tempnn, A.nb, A.mb,
                    mzone, A(m, k), ldam,
                           A(k, n), ldak,
                    zone,  A(m, n), ldam);

                for (m = k+2; m < A.mt; m++)
                {
                    tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                    ldam = BLKLDD(A, m);

                    QUARK_CORE_zgemm_f2(
                        plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempmm, tempnn, A.nb, A.mb,
                        mzone, A(m, k), ldam,
                               A(k, n), ldak,
                        zone,  A(m, n), ldam,
                        /* Dependency on next swapa (gemm need to be done before) */
                        A(k+1, n), A.mb*A.nb, INOUT | GATHERV,
                        /* Dependency on next swapb (gemm need to use panel k before it has to be swaped */
                        fakedep,   1,         INPUT );
                }
            }
        }

        {
            int i, piv[A.nb];
            LAPACKE_zgetrf_work(LAPACK_COL_MAJOR, tempm, tempkn, A0 + tempk * (A.m + 1),  A.m, piv );

            /* Compare ipiv */
            for( i=0; i< min(tempm, tempkn); i++ ) {
                if (piv[i] != (IPIV(k)[i]-tempk) ) {
                    printf("IPIV[%d] = %d - %d (%e - %e)\n", tempk+i, IPIV(k)[i], piv[i]+tempk,
                           A(k, k)[ i * (A.mb+1) ],
                           A0 + (tempk + i) * (A.m + 1) );
                    ok = 0;
                }
            }

            if (ok && (tempk+A.nb < A.n)) {
                LAPACKE_zlaswp_work( LAPACK_COL_MAJOR, A.n - tempk - A.nb,
                                     A0 + (tempk + A.nb) * A.m + tempk, A.m,
                                     1, min(tempm, tempkn), piv, 1);

                CORE_ztrsm(
                    PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaUnit,
                    tempkm, A.n - tempk-A.nb,
                    zone, A0 + tempk * (A.m + 1), A.m,
                          A0 + (tempk + A.nb) * A.m + tempk, A.m);

                if ((tempk+A.mb < A.m)) {
                    CORE_zgemm(
                    PlasmaNoTrans, PlasmaNoTrans,
                    A.m - tempk - A.mb, A.n - tempk - A.nb, A.mb,
                    -1., A0 +  tempk         * A.m + tempk + A.mb, A.m,
                         A0 + (tempk + A.nb) * A.m + tempk,        A.m,
                    1.0, A0 + (tempk + A.nb) * A.m + tempk + A.mb, A.m);
                }
            }
        }
    }

    QUARK_Task_Flag_Set(&task_flags, TASK_PRIORITY, QUARK_TASK_MIN_PRIORITY );
    for (k = 0; k < min(A.mt, A.nt); k++)
    {
        int mintmp;
        tempk  = k * A.mb;
        tempm  = A.m - tempk;
        tempkm = k == A.mt-1 ? tempm : A.mb;
        tempkn = k == A.nt-1 ? A.n - k * A.nb : A.nb;
        mintmp = min(tempkm, tempkn);
        ldak = BLKLDD(A, k);

        /*
         * Apply row interchange behind the panel (work on the panel)
         */
        fakedep = (void*)(intptr_t)k;
        for (n = 0; n < k; n++)
        {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            QUARK_CORE_zlaswp_ontile_f2(
                plasma->quark, &task_flags,
                plasma_desc_submatrix(A, tempk, n*A.nb, tempm, tempnn),
                A(k, n), 1, mintmp, IPIV(k), 1,
                /* Dependency on previous swapb */
                A(k-1,n), A.lm*A.nb, INPUT,
                /* Dependency on all GEMM from previous step */
                fakedep,  1,         INOUT | GATHERV );
        }
    }

    free(A0);
}

