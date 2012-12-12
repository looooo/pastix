/**
 *
 * @file pzshift.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.5.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @precisions normal z -> c d s
 *
 **/

#include <stdlib.h>
#include <sys/types.h>
#include <assert.h>
#include "common.h"
#include "primes.h"
#include "gkkleader.h"

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 *  plasma_zshift Implementation of inplace transposition
 *    based on the GKK algorithm by Gustavson, Karlsson, Kagstrom.
 *    This algorithm shift some cycles to transpose the matrix.
 *
 *******************************************************************************
 *
 * @param[in] plasma
 *         Plasma context
 *
 * @param[in] m
 *         Number of rows of matrix A
 *
 * @param[in] n
 *         Number of columns of matrix A
 *
 * @param[in,out] A
 *         Matrix of size L*m*n
 *
 * @param[in] nprob
 *         Number of parallel and independant problems
 *
 * @param[in] me
 *         Number of rows of the problem
 *
 * @param[in] ne
 *         Number of columns in the problem
 *
 * @param[in] L
 *         Size of chunk to use for transformation
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[in,out] request
 *          Identifies this function call (for exception handling purposes).
 *
 ******************************************************************************/
int plasma_zshift(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A,
                  int nprob, int me, int ne, int L,
                  PLASMA_sequence *sequence, PLASMA_request *request)
{
    int *leaders = NULL;
    int ngrp, thrdbypb, thrdtot, nleaders;

    /* Check Plasma context */
    thrdtot  = PLASMA_SIZE;
    thrdbypb = PLASMA_GRPSIZE;
    ngrp = thrdtot/thrdbypb;

    /* check input */
    if( (nprob * me * ne * L) != (m * n) ) {
        plasma_error(__func__, "problem size does not match matrix size");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if( thrdbypb > thrdtot ) {
        plasma_error(__func__, "number of thread per problem must be less or equal to total number of threads");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if( (thrdtot % thrdbypb) != 0 ) {
        plasma_error(__func__, "number of thread per problem must divide the total number of thread");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /* quick return */
    if( (me < 2) || (ne < 2) || (nprob < 1) ) {
        return PLASMA_SUCCESS;
    }

    /* Get all the cycles leaders and length
     * They are the same for each independent problem (each panel) */
    GKK_getLeaderNbr(me, ne, &nleaders, &leaders);

    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) {
        int *Tp      = NULL;
        int i, ip;
        int owner;
        int *static_leaders;
        int nprob_per_grp = (nprob+ngrp-1)/ngrp;

        Tp = (int *)plasma_shared_alloc(plasma, thrdbypb, PlasmaInteger);
        for (i=0; i<thrdbypb; i++)
            Tp[i] = 0;

        static_leaders = (int *)plasma_shared_alloc(plasma, nleaders*nprob_per_grp*4, PlasmaInteger);
        for (i=0; i<nleaders; i++) {
            static_leaders[i*4  ] = leaders[i*3];
            static_leaders[i*4+1] = leaders[i*3+1];
            static_leaders[i*4+2] = -1;
            static_leaders[i*4+3] = -1;
        }
        for (ip=1; ip<nprob_per_grp; ip++) {
            memcpy(static_leaders + nleaders * 4 * ip,
                   static_leaders, nleaders * 4 * sizeof(int));
        }

        /* loop over leader */
        if (thrdbypb > 1) {
            int *tmp = static_leaders;
            for (ip=0; ip<nprob_per_grp; ip++) {
                for (i=0; i<nleaders; i++) {
                    /* assign this cycle to a thread */
                    owner = minloc(thrdbypb, Tp);

                    /* assign it to owner */
                    Tp[owner] = Tp[owner] + tmp[1] * L;
                    tmp[2] = owner;
                    tmp[3] = ip;
                    tmp += 4;
                }
            }
            /* GKK_BalanceLoad(thrdbypb, Tp, static_leaders, */
            /*                 nleaders, L); */
        }
        else {
            int *tmp = static_leaders;
            for (ip=0; ip<nprob_per_grp; ip++) {
                for (i=0; i<nleaders; i++) {
                    tmp[2] = 0;
                    tmp[3] = ip;
                    tmp += 4;
                }
            }
        }
        nleaders = 4 * nprob_per_grp * nleaders;

        /* shift in parallel */
        plasma_static_call_10(plasma_pzshift,
                              int,                 me,
                              int,                 ne,
                              int,                 L,
                              PLASMA_Complex64_t*, A,
                              int *,               static_leaders,
                              int,                 nleaders,
                              int,                 nprob,
                              int,                 thrdbypb,
                              PLASMA_sequence*,    sequence,
                              PLASMA_request*,     request);

        plasma_shared_free(plasma, Tp);
        plasma_shared_free(plasma, static_leaders);
    }
    /* Dynamic scheduling */
    else {
        nleaders *= 3;
        plasma_dynamic_call_10(plasma_pzshift,
                               int,                 me,
                               int,                 ne,
                               int,                 L,
                               PLASMA_Complex64_t*, A,
                               int *,               leaders,
                               int,                 nleaders,
                               int,                 nprob,
                               int,                 thrdbypb,
                               PLASMA_sequence*,    sequence,
                               PLASMA_request*,     request);
    }

    free(leaders);

    return PLASMA_SUCCESS;
}

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 * plasma_pzshift shifts a batch of cycles in parallel.
 *
 *******************************************************************************
 *
 * @param[in] plasma
 *         Plasma context
 *
 * [in] m
 *         Number of rows of matrix A
 *
 * [in] n
 *         Number of columns of matrix A
 *
 * [in] L
 *         Size of chunk to use for transformation
 *
 * [in,out] A
 *         Matrix of size L*m*n
 *
 * [in] leaders
 *         Array of 4-tuple (cycle leader, cycle length, owner, pb id)
 *
 * [in] nleaders
 *         Number of cycle leaders * 4 * Number of problem per group
 *         Size of the leaders array.
 *
 * [in] nprob
 *         Total Number of parallel and independant problems
 *
 * [in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * [in,out] request
 *          Identifies this function call (for exception handling purposes).
 *
 ******************************************************************************/
void plasma_pzshift(plasma_context_t *plasma) {
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_Complex64_t *A, *Al, *W;
    int     myrank;
    int     i, iprob;
    int     n, m, L, nprob, nleaders, thrdbypb;
    int     locgrp, locrnk, ngrp;
    int    *leaders;

    plasma_unpack_args_10(m, n, L, A, leaders, nleaders, nprob, thrdbypb, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    myrank = PLASMA_RANK;
    locrnk = myrank % thrdbypb;
    locgrp = myrank / thrdbypb;
    ngrp   = PLASMA_SIZE / thrdbypb;

    W = (PLASMA_Complex64_t*)plasma_private_alloc(plasma, L, PlasmaComplexDouble);

    /* shift cycles in parallel. */
    /* each thread shifts the cycles it owns. */
    for(i=0; i<nleaders; i+=4) {
        if( leaders[i+2] == locrnk ) {

            iprob = leaders[i+3] * ngrp + locgrp;

            if ( iprob < nprob ) {
                Al = &(A[iprob*m*n*L]);

                /* cycle #i belongs to this thread, so shift it */
                memcpy(W, &(Al[leaders[i]*L]), L*sizeof(PLASMA_Complex64_t));
                CORE_zshiftw(leaders[i], leaders[i+1], m, n, L, Al, W);
            }
        }
    }

    plasma_private_free(plasma, W);
}


void plasma_pzshift_quark(int m, int n, int L, PLASMA_Complex64_t *A,
                          int *leaders, int nleaders, int nprob, int thrdbypb,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t   *plasma;
    Quark_Task_Flags    task_flags = Quark_Task_Flags_Initializer;
    PLASMA_Complex64_t *Al;
    int     i, iprob, size;
    (void)thrdbypb;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    size = m*n*L;

    for(iprob=0; iprob<nprob; iprob++) {
        Al = &(A[iprob*size]);

        QUARK_Insert_Task(plasma->quark, CORE_foo_quark, &task_flags,
                          sizeof(PLASMA_Complex64_t)*size, Al,  INOUT,
#ifdef TRACE_IPT
                          13, "Foo In shift",   VALUE | TASKLABEL,
                          4, "red",  VALUE | TASKCOLOR,
#endif
                          0);

        /* shift cycles in parallel. */
        for(i=0; i<nleaders; i+=3) {
            //assert( leaders[i+2] != -2 );
            QUARK_CORE_zshift(plasma->quark, &task_flags,
                              leaders[i], m, n, L, Al);
        }

        QUARK_Insert_Task(plasma->quark, CORE_foo_quark, &task_flags,
                          sizeof(PLASMA_Complex64_t)*size, Al,  INOUT,
#ifdef TRACE_IPT
                          14, "Foo Out shift",   VALUE | TASKLABEL,
                          4, "red",  VALUE | TASKCOLOR,
#endif
                          0);
    }
}
