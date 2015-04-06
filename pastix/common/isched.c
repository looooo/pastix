/**
 *
 * @file isched.c
 *
 * Copyright (c) 2008-2014 The University of Bordeaux, IPB, LaBRI, Inria -
 *                         Bordeaux-Sud-Ouest.  All rights reserved.
 *
 * Copyright (c) 2010-2014 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *  PaStiX Internal Thread System routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains basic functions to handle threads for internal schedulings.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "isched.h"

#if defined(HAVE_HWLOC)
#include "isched_hwloc.h"
#elif defined(ARCH_COMPAQ)
#  include <sys/types.h>
#  include <sys/resource.h>
#  include <sys/processor.h>
#  include <sys/sysinfo.h>
#  include <machine/hal_sysinfo.h>
#  define X_INCLUDE_CXML
#elif defined(HAVE_SCHED_SETAFFINITY)
#  include <linux/unistd.h>
#  include <sched.h>
#elif defined(MAC_OS_X)
#  include <mach/mach_init.h>
#  include <mach/thread_policy.h>
/**
 * Expose the hidden kernel interface.
 */
extern kern_return_t thread_policy_set( thread_t               thread,
                                        thread_policy_flavor_t flavor,
                                        thread_policy_t        policy_info,
                                        mach_msg_type_number_t count);
#endif  /* define(HAVE_HWLOC) */

typedef struct isched_barrier_s {
    volatile int    id;
    volatile int    nblocked_thrds;
    pthread_mutex_t synclock;
    pthread_cond_t  synccond;
    int size;
} isched_barrier_t;

typedef struct isched_s {
    isched_barrier_t barrier;

    int              world_size;
    pthread_attr_t   thread_attr;
    int             *rank;
    pthread_t       *tids;

    pthread_mutex_t  statuslock;
    pthread_cond_t   statuscond;
    volatile int     status;

    void           (*pfunc)(void*);
    void            *pargs;
} isched_t;

/***************************************************************************//**
 *  Busy-waiting barrier initialization
 **/
void isched_barrier_init(isched_barrier_t *barrier, int size)
{
    barrier->id = 0;
    barrier->nblocked_thrds = 0;
    barrier->size = size;
    pthread_mutex_init(&(barrier->synclock), NULL);
    pthread_cond_init( &(barrier->synccond), NULL);
}

/***************************************************************************//**
 *  Busy-waiting barrier finalize
 **/
void isched_barrier_finalize(isched_barrier_t *barrier)
{
    pthread_mutex_destroy(&(barrier->synclock));
    pthread_cond_destroy( &(barrier->synccond));
}

/***************************************************************************//**
 *  Non busy-waiting barrier
 **/
void isched_barrier(isched_barrier_t *barrier)
{
    int id;

    pthread_mutex_lock(&(barrier->synclock));
    id = barrier->id;
    barrier->nblocked_thrds++;
    if (barrier->nblocked_thrds == barrier->size) {
        barrier->nblocked_thrds = 0;
        barrier->id++;
        pthread_cond_broadcast(&(barrier->synccond));
    }
    while (id == barrier->id)
        pthread_cond_wait(&(barrier->synccond), &(barrier->synclock));
    pthread_mutex_unlock(&(barrier->synclock));
}

/***************************************************************************//**
 *  Returns core id
 **/
int isched_rank( isched_t *isched )
{
    int rank;
    pthread_t thread_id;

    thread_id = pthread_self();
    for (rank = 0; rank < (isched->barrier).size; rank++)
        if (pthread_equal( isched->tids[rank], thread_id ))
            return rank;
    return -1;
}

/***************************************************************************//**
 *  Main thread control run by each working thread
 **/
void *isched_parallel_section(void *ptr)
{
    isched_t *isched = (isched_t*)(ptr);
    int action;
    int id = isched_rank( isched );
    (void)id;
    assert(id != -1);

    /* Set thread affinity for the worker */
    isched_bind_on_core_index( id );
    isched_barrier( &(isched->barrier) );

    while(1) {
        pthread_mutex_lock( &(isched->statuslock) );
        while ((action = isched->status) == ISCHED_ACT_STAND_BY)
            pthread_cond_wait( &(isched->statuscond), &(isched->statuslock) );
        pthread_mutex_unlock( &(isched->statuslock) );
        isched_barrier( &(isched->barrier) );

        switch (action) {
            case ISCHED_ACT_PARALLEL:
                isched->pfunc( isched->pargs );
                break;
            case ISCHED_ACT_FINALIZE:
                return NULL;
            default:
                fprintf(stderr, "isched_parallel_section: undefined action");
                return NULL;
        }
        isched_barrier(&(isched->barrier) );
    }

    isched_unbind();
    return NULL;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_isched
 *
 *  ischedInit - Initialize the stuctures of the internal thread scheduler.
 *
 *******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use (threads to launch).
 *          If cores = 0, cores = PASTIX_NUM_THREADS if it is set, the
 *          system number of core otherwise.
 *
 * @param[in] coresbind
 *          Array to specify where to bind each thread.
 *          Each thread i is binded to coresbind[hwloc(i)] if hwloc is
 *          provided, or to coresbind[i] otherwise.
 *          If coresbind = NULL, coresbind = PLASMA_AFF_THREADS if it
 *          is set, the identity function otherwise.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS successful exit
 *
 ******************************************************************************/
isched_t *ischedInit(int cores, int *coresbind)
{
    isched_t *isched;
    int status;
    int core;

    /* Create context and insert in the context map */
    MALLOC_INTERN(isched, 1, isched_t);
    if (isched == NULL) {
        fprintf(stderr, "ischedInit: isched allocation failed\n");
        return PASTIX_ERR_OUTOFMEMORY;
    }
    pthread_mutex_init(&(isched->statuslock), NULL);
    pthread_cond_init( &(isched->statuscond), NULL);
    isched->status = ISCHED_ACT_STAND_BY;
    isched->pfunc = NULL;
    isched->pargs = NULL;

    /* Init number of cores and topology */
    isched_init();

    /* Set number of cores */
    if ( cores < 1 ) {
        //isched->world_size = pastix_getenv_int("PASTIX_NUM_THREADS", -1);
        if ( isched->world_size == -1 ) {
            isched->world_size = isched_world_size();
            fprintf(stderr, "ischedInit: Could not find the number of cores: the thread number is set to %d", isched->world_size);
        }
    }
    else
        isched->world_size = cores;

    if (isched->world_size <= 0) {
        fprintf(stderr, "ischedInit: failed to get system size");
        return PASTIX_ERR_INTERNAL;
    }

    /* Initialize barriers */
    isched_barrier_init( &(isched->barrier), isched->world_size );

    /* Initialize default thread attributes */
    status = pthread_attr_init( &(isched->thread_attr) );
    if (status != 0) {
        fprintf(stderr, "ischedInit: pthread_attr_init() failed");
        return status;
    }

    /* Set scope to system */
    status = pthread_attr_setscope( &(isched->thread_attr), PTHREAD_SCOPE_SYSTEM );
    if (status != 0) {
        fprintf(stderr, "ischedInit: pthread_attr_setscope() failed");
        return status;
    }

    /* /\* Set concurrency *\/ */
    /* status = pthread_setconcurrency(isched->world_size); */
    /* if (status != 0) { */
    /*     fprintf(stderr, "ischedInit: pthread_setconcurrency() failed"); */
    /*     return status; */
    /* } */

    /*  Launch threads */
    /* calloc(isched->tids, isched->world_size, sizeof(pthread_t)); */
    /* if (coresbind != NULL) { */
    /*     memcpy(isched->bindings, coresbind, isched->world_size * sizeof(int)); */
    /* } */
    /* else { */
    /*     //plasma_get_affthreads(plasma->thread_bind); */
    /* } */

    MALLOC_INTERN(isched->rank, isched->world_size, int);
    MALLOC_INTERN(isched->tids, isched->world_size, pthread_t);

    /* Assign rank and thread ID for the master */
    isched->rank[0] = 0;
    isched->tids[0] = pthread_self();

    for (core = 1; core < isched->world_size; core++) {
        isched->rank[core] = core;
        pthread_create(
            &isched->tids[core],
            &isched->thread_attr,
             isched_parallel_section,
             (void*)isched);
    }

    return isched;
}

/**
 *****************************************************************************
 *
 * @ingroup pastix_isched
 *
 *  ischedFinalize - Finalize the structures associated to the internal threads
 *  scheduler.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS successful exit
 *
 ******************************************************************************/
int ischedFinalize(isched_t *isched)
{
    int core;
    int status;
    void *exitcodep;

    /* Terminate the dynamic scheduler */
    isched_barrier(&(isched->barrier));

    /* Set termination action */
    pthread_mutex_lock(&isched->statuslock);
    isched->status = ISCHED_ACT_FINALIZE;
    pthread_mutex_unlock(&isched->statuslock);
    pthread_cond_broadcast(&isched->statuscond);

    /* Barrier and clear action */
    isched_barrier(&(isched->barrier));
    isched->status = ISCHED_ACT_STAND_BY;

    // Join threads
    for (core = 1; core < isched->world_size; core++) {
        status = pthread_join(isched->tids[core], &exitcodep);
        if (status != 0) {
            fprintf(stderr, "ischedFinalize: pthread_join() failed");
            return status;
        }
    }
    isched_barrier_finalize(&(isched->barrier));

    /* Unbind main thread */
    isched_unbind();

    /* Destroy thread attributes */
    status = pthread_attr_destroy(&isched->thread_attr);
    if (status != 0)
        fprintf(stderr, "ischedFinalize: pthread_attr_destroy() failed");

    /* Destroy topology */
    isched_finalize();

    memFree_null(isched->rank);
    memFree_null(isched->tids);

    memFree_null(isched);
    return PASTIX_SUCCESS;
}
