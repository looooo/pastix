/**
 *
 * @file control.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#if defined( _WIN32 ) || defined( _WIN64 )
#include "plasmawinthread.h"
#else
#include <pthread.h>
#endif
#include "common.h"
#include "auxiliary.h"
#include "allocate.h"

/***************************************************************************//**
 *  Busy-waiting barrier initialization
 **/
void plasma_barrier_init(plasma_context_t *plasma)
{
#ifdef BUSY_WAITING
    int core;

    for (core = 0; core < CONTEXT_THREADS_MAX; core++) {
        plasma->barrier_in[core] = 0;
        plasma->barrier_out[core] = 0;
    }
#else
    plasma->barrier_id = 0;
    plasma->barrier_nblocked_thrds = 0;
    pthread_mutex_init(&(plasma->barrier_synclock), NULL);
    pthread_cond_init( &(plasma->barrier_synccond), NULL);
#endif
}

/***************************************************************************//**
 *  Busy-waiting barrier finalize
 **/
void plasma_barrier_finalize(plasma_context_t *plasma)
{
#ifndef BUSY_WAITING
    pthread_mutex_destroy(&(plasma->barrier_synclock));
    pthread_cond_destroy( &(plasma->barrier_synccond));
#endif
}

/***************************************************************************//**
 *  Busy-waiting barrier
 **/
void plasma_barrier(plasma_context_t *plasma)
{
#ifdef BUSY_WAITING
    int core;

    if (PLASMA_RANK == 0) {
        for (core = 1; core < PLASMA_SIZE; core++)
            while (plasma->barrier_in[core] == 0);

        for (core = 1; core < PLASMA_SIZE; core++)
            plasma->barrier_in[core] = 0;

        for (core = 1; core < PLASMA_SIZE; core++)
            plasma->barrier_out[core] = 1;
    }
    else
    {
        plasma->barrier_in[PLASMA_RANK] = 1;
        while (plasma->barrier_out[PLASMA_RANK] == 0);
        plasma->barrier_out[PLASMA_RANK] = 0;
    }
#else
    int id;

    pthread_mutex_lock(&(plasma->barrier_synclock));
    id = plasma->barrier_id;
    plasma->barrier_nblocked_thrds++;
    if (plasma->barrier_nblocked_thrds == PLASMA_SIZE) {
        plasma->barrier_nblocked_thrds = 0;
        plasma->barrier_id++;
        pthread_cond_broadcast(&(plasma->barrier_synccond));
    }
    while (id == plasma->barrier_id)
        pthread_cond_wait(&(plasma->barrier_synccond), &(plasma->barrier_synclock));
    pthread_mutex_unlock(&(plasma->barrier_synclock));
#endif
}

/***************************************************************************//**
 *  Main thread control
 **/
void *plasma_parallel_section(void *plasma_ptr)
{
    plasma_context_t *plasma = (plasma_context_t*)(plasma_ptr);
    PLASMA_enum action;

    /* Set thread affinity for the worker */
    plasma_setaffinity(plasma->thread_bind[plasma_rank(plasma)]);

    plasma_barrier(plasma);
    while(1) {
        pthread_mutex_lock(&plasma->action_mutex);
        while ((action = plasma->action) == PLASMA_ACT_STAND_BY)
            pthread_cond_wait(&plasma->action_condt, &plasma->action_mutex);
        pthread_mutex_unlock(&plasma->action_mutex);
        plasma_barrier(plasma);

        switch (action) {
            case PLASMA_ACT_PARALLEL:
                plasma->parallel_func_ptr(plasma);
                break;
            case PLASMA_ACT_DYNAMIC:
                QUARK_Worker_Loop(plasma->quark, plasma_rank(plasma));
                break;
            case PLASMA_ACT_FINALIZE:
                return NULL;
            default:
                plasma_fatal_error("plasma_parallel_section", "undefined action");
                return NULL;
        }
        plasma_barrier(plasma);
    }

    plasma_unsetaffinity();
    return NULL;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Init - Initialize PLASMA.
 *
 *******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use (threads to launch).
 *          If cores = 0, cores = PLASMA_NUM_THREADS if it is set, the
 *          system number of core otherwise.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Init(int cores)
{
    return PLASMA_Init_Affinity(cores, NULL);
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Init_Affinity - Initialize PLASMA.
 *
 *******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use (threads to launch).
 *          If cores = 0, cores = PLASMA_NUM_THREADS if it is set, the
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
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Init_Affinity(int cores, int *coresbind)
{
    plasma_context_t *plasma;
    int status;
    int core;

    /* Create context and insert in the context map */
    plasma = plasma_context_create();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Init", "plasma_context_create() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    status = plasma_context_insert(plasma, pthread_self());
    if (status != PLASMA_SUCCESS) {
        plasma_fatal_error("PLASMA_Init", "plasma_context_insert() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    /* Init number of cores and topology */
    plasma_topology_init();

    /* Set number of cores */
    if ( cores < 1 ) {
        plasma->world_size = plasma_get_numthreads();
        if ( plasma->world_size == -1 ) {
            plasma->world_size = 1;
            plasma_warning("PLASMA_Init", "Could not find the number of cores: the thread number is set to 1");
        }
    }
    else
      plasma->world_size = cores;

    if (plasma->world_size <= 0) {
        plasma_fatal_error("PLASMA_Init", "failed to get system size");
        return PLASMA_ERR_NOT_FOUND;
    }
    /* Check if not more cores than the hard limit */
    if (plasma->world_size > CONTEXT_THREADS_MAX) {
        plasma_fatal_error("PLASMA_Init", "not supporting so many cores");
        return PLASMA_ERR_INTERNAL_LIMIT;
    }

    /* Get the size of each NUMA node */
    plasma->group_size = plasma_get_numthreads_numa();
    while ( ((plasma->world_size)%(plasma->group_size)) != 0 ) 
        (plasma->group_size)--;

    /* Initialize barrier */
    plasma_barrier_init(plasma);

    /* Initialize default thread attributes */
    status = pthread_attr_init(&plasma->thread_attr);
    if (status != 0) {
        plasma_fatal_error("PLASMA_Init", "pthread_attr_init() failed");
        return status;
    }
    /* Set scope to system */
    status = pthread_attr_setscope(&plasma->thread_attr, PTHREAD_SCOPE_SYSTEM);
    if (status != 0) {
        plasma_fatal_error("PLASMA_Init", "pthread_attr_setscope() failed");
        return status;
    }
    /* Set concurrency */
    status = pthread_setconcurrency(plasma->world_size);
    if (status != 0) {
        plasma_fatal_error("PLASMA_Init", "pthread_setconcurrency() failed");
        return status;
    }
    /*  Launch threads */
    memset(plasma->thread_id,   0, CONTEXT_THREADS_MAX*sizeof(pthread_t));
    if (coresbind != NULL) {
        memcpy(plasma->thread_bind, coresbind, plasma->world_size*sizeof(int));
    }
    else {
        plasma_get_affthreads(plasma->thread_bind);
    }
    /* Assign rank and thread ID for the master */
    plasma->thread_rank[0] = 0;
    plasma->thread_id[0] = pthread_self();

    for (core = 1; core < plasma->world_size; core++) {
        plasma->thread_rank[core] = core;
        pthread_create(
            &plasma->thread_id[core],
            &plasma->thread_attr,
             plasma_parallel_section,
             (void*)plasma);
    }
    /* Set thread affinity for the master */
    plasma_setaffinity(plasma->thread_bind[0]);

    /* Initialize the dynamic scheduler */
    plasma->quark =  QUARK_Setup(plasma->world_size);

    plasma_barrier(plasma);
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Finalize - Finalize PLASMA.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Finalize()
{
    int core;
    int status;
    void *exitcodep;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Finalize()", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Terminate the dynamic scheduler */
    plasma_dynamic_sync();

    /* Free quark structures */
    QUARK_Free(plasma->quark);

    /* Set termination action */
    pthread_mutex_lock(&plasma->action_mutex);
    plasma->action = PLASMA_ACT_FINALIZE;
    pthread_mutex_unlock(&plasma->action_mutex);
    pthread_cond_broadcast(&plasma->action_condt);

    /* Barrier and clear action */
    plasma_barrier(plasma);
    plasma->action = PLASMA_ACT_STAND_BY;

    // Join threads
    for (core = 1; core < plasma->world_size; core++) {
        status = pthread_join(plasma->thread_id[core], &exitcodep);
        if (status != 0) {
            plasma_fatal_error("PLASMA_Finalize", "pthread_join() failed");
            return status;
        }
    }
    plasma_barrier_finalize(plasma);

    /* Unbind main thread */
    plasma_unsetaffinity();

    /* Destroy thread attributes */
    status = pthread_attr_destroy(&plasma->thread_attr);
    if (status != 0)
        plasma_fatal_error("PLASMA_Finalize", "pthread_attr_destroy() failed");

    /* Destroy topology */
    plasma_topology_finalize();

    status = plasma_context_remove(plasma, pthread_self());
    if (status != PLASMA_SUCCESS) {
        plasma_fatal_error("PLASMA_Finalize", "plasma_context_remove() failed");
        return status;
    }

    /* Restore the concurency */
    /* actually it's really bad, we shoulde set the concurrency only
     * if it's not already done and restore it only we had change it */
    pthread_setconcurrency( 0 );

    return PLASMA_SUCCESS;
}
