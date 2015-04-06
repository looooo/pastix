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

    pthread_attr_t   thread_attr;
    pthread_mutex_t  statuslock;
    pthread_cond_t   statuscond;
    volatile int     status;
    pthread_t       *tids;
    void           (*pfunc)(void*);
    void            *pargs;
} isched_t;

/***************************************************************************//**
 *  Busy-waiting barrier initialization
 **/
void isched_barrier_init(isched_barrier_t *barrier)
{
    barrier->id = 0;
    barrier->nblocked_thrds = 0;
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

/* int ischedInit(int cores, int *coresbind) */
/* { */
/*     isched_t *isched; */
/*     int status; */
/*     int core; */

/*     /\* Create context and insert in the context map *\/ */
/*     isched = (isched_t*)malloc(sizeof(isched_t)); */
/*     if (isched == NULL) { */
/*         fprintf(stderr, "ischedInit: isched allocation failed\n"); */
/*         return PASTIX_ERR_OUT_OF_RESOURCES; */
/*     } */
/*     status = plasma_context_insert(plasma, pthread_self()); */
/*     if (status != PLASMA_SUCCESS) { */
/*         plasma_fatal_error("PLASMA_Init", "plasma_context_insert() failed"); */
/*         return PLASMA_ERR_OUT_OF_RESOURCES; */
/*     } */
/*     /\* Init number of cores and topology *\/ */
/*     plasma_topology_init(); */

/*     /\* Set number of cores *\/ */
/*     if ( cores < 1 ) { */
/*         plasma->world_size = plasma_get_numthreads(); */
/*         if ( plasma->world_size == -1 ) { */
/*             plasma->world_size = 1; */
/*             plasma_warning("PLASMA_Init", "Could not find the number of cores: the thread number is set to 1"); */
/*         } */
/*     } */
/*     else */
/*       plasma->world_size = cores; */

/*     if (plasma->world_size <= 0) { */
/*         plasma_fatal_error("PLASMA_Init", "failed to get system size"); */
/*         return PLASMA_ERR_NOT_FOUND; */
/*     } */
/*     /\* Check if not more cores than the hard limit *\/ */
/*     if (plasma->world_size > CONTEXT_THREADS_MAX) { */
/*         plasma_fatal_error("PLASMA_Init", "not supporting so many cores"); */
/*         return PLASMA_ERR_INTERNAL_LIMIT; */
/*     } */

/*     /\* Get the size of each NUMA node *\/ */
/*     plasma->group_size = plasma_get_numthreads_numa(); */
/*     while ( ((plasma->world_size)%(plasma->group_size)) != 0 ) */
/*         (plasma->group_size)--; */

/*     /\* Initialize barriers *\/ */
/*     isched_barrier_init(&(isched->barrier)); */

/*     /\* Initialize default thread attributes *\/ */
/*     status = pthread_attr_init(&plasma->thread_attr); */
/*     if (status != 0) { */
/*         plasma_fatal_error("PLASMA_Init", "pthread_attr_init() failed"); */
/*         return status; */
/*     } */
/*     /\* Set scope to system *\/ */
/*     status = pthread_attr_setscope(&plasma->thread_attr, PTHREAD_SCOPE_SYSTEM); */
/*     if (status != 0) { */
/*         plasma_fatal_error("PLASMA_Init", "pthread_attr_setscope() failed"); */
/*         return status; */
/*     } */
/*     /\* Set concurrency *\/ */
/*     status = pthread_setconcurrency(plasma->world_size); */
/*     if (status != 0) { */
/*         plasma_fatal_error("PLASMA_Init", "pthread_setconcurrency() failed"); */
/*         return status; */
/*     } */
/*     /\*  Launch threads *\/ */
/*     memset(plasma->thread_id,   0, CONTEXT_THREADS_MAX*sizeof(pthread_t)); */
/*     if (coresbind != NULL) { */
/*         memcpy(plasma->thread_bind, coresbind, plasma->world_size*sizeof(int)); */
/*     } */
/*     else { */
/*         plasma_get_affthreads(plasma->thread_bind); */
/*     } */
/*     /\* Assign rank and thread ID for the master *\/ */
/*     plasma->thread_rank[0] = 0; */
/*     plasma->thread_id[0] = pthread_self(); */

/*     for (core = 1; core < plasma->world_size; core++) { */
/*         plasma->thread_rank[core] = core; */
/*         pthread_create( */
/*             &plasma->thread_id[core], */
/*             &plasma->thread_attr, */
/*              plasma_parallel_section, */
/*              (void*)plasma); */
/*     } */

/*     /\* Ensure BLAS are sequential and set thread affinity for the master *\/ */
/* #if defined(PLASMA_WITH_MKL) */
/* #if defined(__ICC) || defined(__INTEL_COMPILER) */
/*     kmp_set_defaults("KMP_AFFINITY=disabled"); */
/* #endif */
/* #endif */

/*     /\* Initialize the dynamic scheduler *\/ */
/*     plasma->quark =  QUARK_Setup(plasma->world_size); */
/*     plasma_barrier(plasma); */

/*     plasma_setlapack_sequential(plasma); */

/*     return PLASMA_SUCCESS; */
/* } */
