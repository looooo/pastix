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
 * @version 2.4.6
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

    GKK_getLeaderNbr(me, ne, &nleaders, &leaders);
    nleaders *= 3;

    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) {
        int *Tp      = NULL;
        int i, ipb;
        int owner;

        Tp = (int *)plasma_shared_alloc(plasma, thrdtot, PlasmaInteger);
        for (i=0; i<thrdtot; i++)
            Tp[i] = 0;

        ipb = 0;
        
        /* First part with coarse parallelism */
        if (nprob > ngrp) {
            ipb = (nprob / ngrp)*ngrp;
        
            /* loop over leader */
            if (thrdbypb > 1) {
                for (i=0; i<nleaders; i+=3) {
                    /* assign this cycle to a thread */
                    owner = minloc(thrdbypb, Tp);
                
                    /* assign it to owner */
                    Tp[owner] = Tp[owner] + leaders[i+1] * L;
                    leaders[i+2] = owner;
                }
            
                GKK_BalanceLoad(thrdbypb, Tp, leaders, nleaders, L);
            }
            else {
                for (i=0; i<nleaders; i+=3) {
                    Tp[0] = Tp[0] + leaders[i+1] * L;
                    leaders[i+2] = 0;
                }
            }

            /* shift in parallel */
            for (i=0; i< (nprob/ngrp); i++) {
                plasma_static_call_9(plasma_pzshift,
                                     int,                 me,
                                     int,                 ne,
                                     int,                 L,
                                     PLASMA_Complex64_t*, &(A[i*ngrp*me*ne*L]),
                                     int *,               leaders,
                                     int,                 nleaders,
                                     int,                 thrdbypb,
                                     PLASMA_sequence*,    sequence,
                                     PLASMA_request*,     request);
            }
        }
    
        /* Second part with fine parallelism */
        if (ipb < nprob) {
            for (i=0; i<thrdtot; i++)
                Tp[i] = 0;
        
            if (thrdtot > 1) {
                /* loop over leader */
                for (i=0; i<nleaders; i+=3) {
                    /* assign this cycle to a thread */
                    owner = minloc(thrdtot, Tp);
                
                    /* assign it to owner */
                    Tp[owner] = Tp[owner] + leaders[i+1] * L;
                    leaders[i+2] = owner;
                }
                GKK_BalanceLoad(thrdtot, Tp, leaders, nleaders, L);
            }
            else {
                for (i=0; i<nleaders; i+=3) {
                    Tp[0] = Tp[0] + leaders[i+1] * L;
                    leaders[i+2] = 0;
                }
            }
        
            /* shift in parallel */
            for (i=ipb; i<nprob; i++) {
                plasma_static_call_9(plasma_pzshift,
                                     int,                 me,
                                     int,                 ne,
                                     int,                 L,
                                     PLASMA_Complex64_t*, &(A[i*me*ne*L]),
                                     int *,               leaders,
                                     int,                 nleaders,
                                     int,                 thrdtot,
                                     PLASMA_sequence*,    sequence,
                                     PLASMA_request*,     request);
            }
        }

        plasma_shared_free(plasma, Tp);
    }
    /* Dynamic scheduling */
    else {
        plasma_dynamic_call_9(plasma_pzshift,
                              int,                 me,
                              int,                 ne,
                              int,                 L,
                              PLASMA_Complex64_t*, A,
                              int *,               leaders,
                              int,                 nleaders,
                              int,                 nprob,
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
 * [in,out] A
 *         Matrix of size L*m*n
 *
 * [in] nprob
 *         Number of parallel and independant problems
 *
 * [in] me
 *         Number of rows of the problem
 *
 * [in] ne
 *         Number of columns in the problem
 *
 * [in] L
 *         Size of chunk to use for transformation
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
    int     locrnk, myrank;
    int     i, x, snix, cl, iprob;
    int     n, m, L, nleaders, thrdbypb;
    int    *leaders;
    int64_t s, q;
    
    plasma_unpack_args_9(m, n, L, A, leaders, nleaders, thrdbypb, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    myrank   = PLASMA_RANK;
    locrnk   = myrank % thrdbypb;
    iprob    = myrank / thrdbypb;

    q  = m * n - 1;
    Al = &(A[iprob*m*n*L]);
    
    W = (PLASMA_Complex64_t*)plasma_private_alloc(plasma, L, PlasmaComplexDouble);

    /* shift cycles in parallel. */
    /* each thread shifts the cycles it owns. */
    for(i=0; i<nleaders; i+=3) {
        if( leaders[i+2] == locrnk ) {
            /* cycle #i belongs to this thread, so shift it */
            memcpy(W, &(Al[leaders[i]*L]), L*sizeof(PLASMA_Complex64_t));
            CORE_zshiftw(leaders[i], leaders[i+1], m, n, L, Al, W);
        }
        else if( leaders[i+2] == -2 ) {
            /* cycle #i has been split, so shift in parallel */
            x  = leaders[i+1] / thrdbypb;
            cl = x;
            if( locrnk == 0 ) {
                cl = leaders[i+1] - x * (thrdbypb - 1);
            }
            s    = leaders[i];
            snix = (s * modpow(n, locrnk*x, m * n - 1)) % q;
            
            /* copy the block at s*n^(thid*x) (snix) */
            memcpy(W, &(Al[snix*L]), L*sizeof(PLASMA_Complex64_t));

            /* wait for peers to finish copy their block. */
            plasma_barrier(plasma);

            /* shift the linear array. */
            if( cl > 0 ) {
                CORE_zshiftw(snix, cl, m, n, L, Al, W);
            }
        }
    }

    plasma_private_free(plasma, W);
}


void plasma_pzshift_quark(int m, int n, int L, PLASMA_Complex64_t *A, 
                          int *leaders, int nleaders, int nprob,
                          PLASMA_sequence *sequence, PLASMA_request *request) 
{
    plasma_context_t   *plasma;
    Quark_Task_Flags    task_flags = Quark_Task_Flags_Initializer;
    PLASMA_Complex64_t *Al;
    int     i, iprob, size;
    
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


