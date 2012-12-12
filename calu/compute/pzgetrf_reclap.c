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
 * @version 2.5.0
 * @author Mathieu Faverge
 * @author Hatem Ltaief
 * @date 2009-11-15
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

void CORE_zgetrf_reclap_init(void);

#define PARALLEL_KERNEL
#define LAPACK_LAYOUT
#ifdef LAPACK_LAYOUT
#undef BLKLDD
#define BLKLDD(A, k) (A).lm
#define A(m,n) (&((PLASMA_Complex64_t*)(A.mat))[(int64_t)(A.lm)*(int64_t)(A.nb)*(int64_t)(n)+(int64_t)(A.mb)*(int64_t)(m)])
#else
#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#endif

#define IPIV(k) &(IPIV[(int64_t)A.mb*(int64_t)(k)])

#define plasma_pzgetrf_reclap_rl_quark plasma_pzgetrf_reclap_quark

/***************************************************************************//**
 *  Parallel tile LU factorization - dynamic scheduling - Right looking
 **/
void plasma_pzgetrf_reclap_rl_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request)
{
    int k, m, n, minmnt;
    plasma_context_t *plasma;
    int tempkm, tempkn, tempmm, tempnn;
    int tempm, tempk;
    int ldak, ldam;
    Quark_Task_Flags task_flagsP = Quark_Task_Flags_Initializer;
    Quark_Task_Flags task_flagsU = Quark_Task_Flags_Initializer;

    PLASMA_Complex64_t zone  = (PLASMA_Complex64_t)1.0;
    PLASMA_Complex64_t mzone = (PLASMA_Complex64_t)-1.0;

    void * fakedep;
     /* How many threads per panel? Probably needs to be adjusted during factorization. */
    int panel_thread_count;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flagsP, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    QUARK_Task_Flag_Set(&task_flagsU, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    /* We start at PLASMA_SIZE-1, to keep the first thread adding task to the queue */
    /* kernel doesn't accept more than 48 cores */
    panel_thread_count = min(PLASMA_SIZE-1, 48);

    QUARK_Task_Flag_Set(&task_flagsP, TASK_THREAD_COUNT, panel_thread_count );

    CORE_zgetrf_reclap_init();

    minmnt = min(A.mt, A.nt);
    for (k = 0; k < minmnt; k++)
    {
        tempk  = k * A.mb;
        tempm  = A.m - k * A.mb;
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        ldak = BLKLDD(A, k);

        QUARK_Task_Flag_Set(&task_flagsU, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - k );
#ifdef PARALLEL_KERNEL
        while ( (panel_thread_count * 4 * A.mb) > tempm ) {
            panel_thread_count--;
            QUARK_Task_Flag_Set(&task_flagsP, TASK_THREAD_COUNT, panel_thread_count );
        }
          
        if ( panel_thread_count > 1 ) {
            QUARK_Task_Flag_Set(&task_flagsP, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - k );
            
            QUARK_CORE_zgetrf_reclap(
                plasma->quark, &task_flagsP,
                tempm, tempkn, A.nb,
                A(k, k), ldak, IPIV(k),
                sequence, request, 1, tempk,
                panel_thread_count );
        }
        else {
            QUARK_CORE_zgetrf(
                plasma->quark, &task_flagsU,
                tempm, tempkn, A.mb,
                A(k, k), ldak, IPIV(k),
                sequence, request, 1, tempk );
        }
#else
        QUARK_CORE_zgetrf(
                plasma->quark, &task_flagsU,
                tempm, tempkn, A.mb,
                A(k, k), ldak, IPIV(k),
                sequence, request, 1, tempk );
#endif

        fakedep = (void *)(intptr_t)(k+1);
        for (n = k+1; n < A.nt; n++)
        {

            QUARK_Task_Flag_Set(&task_flagsU, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - n );
            /*
             * Apply row interchange after the panel (work on the panel)
             */
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            QUARK_CORE_zlaswp(
                plasma->quark, &task_flagsU,
                tempnn, A(k, n), A.lm, 1, tempkm, IPIV(k), 1);

            QUARK_CORE_ztrsm(
                plasma->quark, &task_flagsU,
                PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaUnit,
                tempkm, tempnn, A.mb,
                zone, A(k, k), ldak,
                      A(k, n), ldak);

            m = k+1;
            if ( m < A.mt ) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                
                QUARK_CORE_zgemm2(
                    plasma->quark, &task_flagsU,
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
                        plasma->quark, &task_flagsU,
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
    }

    for (k = 0; k < minmnt; k++)
    {
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        tempk  = min(tempkn, tempkm);
        ldak = BLKLDD(A, k);

        /*
         * Apply row interchange behind the panel (work on the panel)
         */
        QUARK_Task_Flag_Set(&task_flagsU, TASK_PRIORITY, QUARK_TASK_MIN_PRIORITY );
        fakedep = (void*)(intptr_t)k;
        for (n = 0; n < k; n++)
        {
            QUARK_CORE_zlaswp_f2(
                plasma->quark, &task_flagsU,
                A.nb, A(k, n), A.lm, 1, tempk, IPIV(k), 1,
                /* Dependency on previous swapb */
                A(k-1,n), A.lm*A.nb, INPUT,
                /* Dependency on all GEMM from previous step */
                fakedep,  1,         INOUT | GATHERV );
        }
    }
}

/***************************************************************************//**
 *  Parallel tile LU factorization - dynamic scheduling - Left looking
 **/
void plasma_pzgetrf_reclap_ll_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request)
{
    int k, m, n;
    plasma_context_t *plasma;
    int tempkm, tempkn, tempmm, tempnn;
    int tempm;
    int ldak, ldam;
    Quark_Task_Flags task_flagsP = Quark_Task_Flags_Initializer;
    Quark_Task_Flags task_flagsU = Quark_Task_Flags_Initializer;

    PLASMA_Complex64_t zone  = (PLASMA_Complex64_t)1.0;
    PLASMA_Complex64_t mzone = (PLASMA_Complex64_t)-1.0;

    void * fakedep;
    int panel_thread_count; /* How many threads per panel? Probably needs to be adjusted during factorization. */

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flagsP, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    QUARK_Task_Flag_Set(&task_flagsU, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    /* We start at PLASMA_SIZE-1, to keep the first thread adding task to the queue */
    panel_thread_count = min(PLASMA_SIZE-1, 48); /* kernel doesn't accept more than 48 cores */
    QUARK_Task_Flag_Set(&task_flagsP, TASK_THREAD_COUNT, panel_thread_count );

    CORE_zgetrf_reclap_init();

    fakedep = (void*)(intptr_t)1;
    for (n = 0; n < A.nt; n++)
    {
        tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;

        QUARK_Task_Flag_Set(&task_flagsU, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - n );
        QUARK_Task_Flag_Set(&task_flagsP, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - n );

        for (k = 0; k < min(A.mt, n); k++)
        {
            tempm  = A.m - k * A.mb;
            tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
            ldak = BLKLDD(A, k);
            
            QUARK_CORE_zlaswp(
                plasma->quark, &task_flagsU,
                tempnn, A(k, n), A.lm, 1, tempkm, IPIV(k), 1);

            QUARK_CORE_ztrsm(
                plasma->quark, &task_flagsU,
                PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaUnit,
                tempkm, tempnn, A.mb,
                zone, A(k, k), ldak,
                      A(k, n), ldak);

            if (k < A.mt-1) {
                m = k+1;
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);

                QUARK_CORE_zgemm2(
                    plasma->quark, &task_flagsU,
                    PlasmaNoTrans, PlasmaNoTrans,
                    tempmm, tempnn, A.nb, A.mb,
                    mzone, A(m, k), ldam,
                           A(k, n), ldak,
                    zone,  A(m, n), ldam);

                fakedep = (void*)(intptr_t)k;
                for (m = k+2; m < A.mt; m++)
                {
                    tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                    ldam = BLKLDD(A, m);
                    
                    QUARK_CORE_zgemm_f2(
                        plasma->quark, &task_flagsU,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempmm, tempnn, A.nb, A.mb,
                        mzone, A(m, k), ldam,
                        A(k, n), ldak,
                        zone,  A(m, n), ldam,
                        /* Dependency on next swapa or getrf (gemm need to be done before) */
                        A(k+1, n), A.mb*A.nb, INOUT | GATHERV,
                        /* Dependency on next swapb (gemm need to use panel k before it has to be swaped */
                        fakedep,   1,         INPUT );
                }
            }
        }

        k = n;
        if ( n < A.mt ) {
            tempm  = A.m - k * A.mb;
            tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
            ldak = BLKLDD(A, k);

#ifdef PARALLEL_KERNEL
            while ( (panel_thread_count * 4 * A.mb + 1) > tempm ) {
                panel_thread_count = panel_thread_count >> 1;
                QUARK_Task_Flag_Set(&task_flagsP, TASK_THREAD_COUNT, panel_thread_count );
            }

            if ( panel_thread_count > 1 ) {
                QUARK_CORE_zgetrf_reclap(
                    plasma->quark, &task_flagsP,
                    tempm, tempkn, A.mb,
                    A(k, k), ldak, IPIV(k),
                    sequence, request, 1, A.mb*k,
                    panel_thread_count );
            } else {
                QUARK_CORE_zgetrf(
                    plasma->quark, &task_flagsU,
                    tempm, tempkn, A.mb,
                    A(k, k), ldak, IPIV(k),
                    sequence, request, 1, A.mb*k );
            }
#else
            QUARK_CORE_zgetrf(
                plasma->quark, &task_flagsU,
                tempm, tempkn, A.mb,
                A(k, k), ldak, IPIV(k),
                sequence, request, 1, A.mb*k );
#endif
        }
    }


    QUARK_Task_Flag_Set(&task_flagsU, TASK_PRIORITY, QUARK_TASK_MIN_PRIORITY );
    for (k = 0; k < min(A.mt, A.nt); k++)
    {
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        ldak = BLKLDD(A, k);
        
        fakedep = (void*)(intptr_t)k;
        for (n = 0; n < k; n++)
        {
            /*
             * Apply row interchange behind the panel (work on the panel)
             */
            QUARK_CORE_zlaswp_f2(
                plasma->quark, &task_flagsU,
                A.nb, A(k, n), ldak, 1, tempkm, IPIV(k), 1,
                /* Dependency on previous swapb */
                A(k-1, n), A.lm*A.nb, INPUT,
                /* Dependency on all GEMM from previous step */
                fakedep,  1,         INOUT | GATHERV );
        }
    }
}
