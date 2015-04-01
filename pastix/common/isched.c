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
    pthread_mutex_t  statuslock;
    pthread_cond_t   statuscond;
    volatile int     status;
    pthread_t       *tids;
    void           (*pfunc)(void*);
    void            *pargs;
} isched_t;

int isched_setaffinity(int cpu)
{
#if defined(HAVE_HWLOC)
    {
        cpu = isched_hwloc_bind_on_core_index(cpu);
        if(cpu == -1 ) {
            fprintf(stderr, "Core binding on node %i failed\n", cpu);
            return -1;
        }
    }
#else /* We bind thread ourself in funtion of the targetted architecture */

#if defined(HAVE_SCHED_SETAFFINITY)
    {
        cpu_set_t mask;
        CPU_ZERO(&mask);
        CPU_SET(cpu, &mask);

#if defined(HAVE_OLD_SCHED_SETAFFINITY)
        if(sched_setaffinity(0,&mask) < 0)
#else /* HAVE_OLD_SCHED_SETAFFINITY */
        if(sched_setaffinity(0,sizeof(mask),&mask) < 0)
#endif /* HAVE_OLD_SCHED_SETAFFINITY */
        {
            return -1;
        }
    }
#elif defined(ARCH_PPC)
    {
        tid_t self_ktid = thread_self();
        bindprocessor(BINDTHREAD, self_ktid, cpu*2);
    }
#elif defined(ARCH_COMPAQ)
    {
        bind_to_cpu_id(getpid(), cpu, 0);
    }
#elif defined(MAC_OS_X)
    {
        thread_affinity_policy_data_t ap;
        int                           ret;

        ap.affinity_tag = 1; /* non-null affinity tag */
        ret = thread_policy_set(
                                mach_thread_self(),
                                THREAD_AFFINITY_POLICY,
                                (integer_t*) &ap,
                                THREAD_AFFINITY_POLICY_COUNT
                                );
        if(ret != 0) {
            return -1;
        }
    }
#endif /* Architectures */

#endif /* WITH_HWLOC     */

    return cpu;
}

int isched_unsetaffinity()
{
#if defined(HAVE_HWLOC)
    {
        return isched_hwloc_unbind();
    }
#else /* We bind thread ourself in funtion of the targetted architecture */
#ifndef PLASMA_AFFINITY_DISABLE
#if (defined PLASMA_OS_LINUX) || (defined PLASMA_OS_FREEBSD)
    {
        int i;
#if (defined PLASMA_OS_LINUX)
        cpu_set_t set;
#elif (defined PLASMA_OS_FREEBSD)
        cpuset_t set;
#endif
        CPU_ZERO( &set );

        for(i=0; i<sys_corenbr; i++)
            CPU_SET( i, &set );

#if (defined HAVE_OLD_SCHED_SETAFFINITY)
        if( sched_setaffinity( 0, &set ) < 0 )
#else /* HAVE_OLD_SCHED_SETAFFINITY */
#if (defined PLASMA_OS_LINUX)
        if( sched_setaffinity( 0, sizeof(set), &set) < 0 )
#elif (defined PLASMA_OS_FREEBSD)
        if( cpuset_setaffinity(CPU_LEVEL_WHICH, CPU_WHICH_PID, 0, sizeof(set), &set) < 0 )
#endif
#endif /* HAVE_OLD_SCHED_SETAFFINITY */
            {
                plasma_warning("plasma_unsetaffinity", "Could not unbind thread");
                return PLASMA_ERR_UNEXPECTED;
            }

        return PLASMA_SUCCESS;
    }
#elif (defined PLASMA_OS_MACOS)
    {
        /* TODO: check how to unbind the main thread if necessary for OpenMP */
        /* thread_affinity_policy_data_t ap; */
        /* int                           ret; */

        /* ap.affinity_tag = 1; /\* non-null affinity tag *\/ */
        /* ret = thread_policy_set( mach_thread_self(), */
        /*                          THREAD_AFFINITY_POLICY, */
        /*                          (integer_t*) &ap, */
        /*                          THREAD_AFFINITY_POLICY_COUNT */
        /*     ); */
        /* if(ret != 0) { */
        /*     plasma_warning("plasma_unsetaffinity", "Could not unbind thread"); */
        /*     return PLASMA_ERR_UNEXPECTED; */
        /* } */

        return PLASMA_SUCCESS;
    }
#elif (defined PLASMA_OS_WINDOWS)
    {
        int i;
        DWORD mask = 0;

        for(i=0; i<sys_corenbr; i++)
            mask |= 1 << i;

        if( SetThreadAffinityMask(GetCurrentThread(), mask) == 0) {
            plasma_warning("plasma_unsetaffinity", "Could not unbind thread");
            return PLASMA_ERR_UNEXPECTED;
        }
        return PLASMA_SUCCESS;
    }
#elif (defined PLASMA_OS_AIX)
    {
        /* TODO: check how to unbind the main thread if necessary for OpenMP */
        /* tid_t self_ktid = thread_self (); */
        /* bindprocessor(BINDTHREAD, self_ktid, rank); */
        return PLASMA_SUCCESS;
    }
#else
    return PLASMA_ERR_NOT_SUPPORTED;
#endif
#endif /* PLASMA_AFFINITY_DISABLE */
#endif /* PLASMA_HWLOC */
}

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
    //isched_setaffinity( pastix->bindtab[isched_rank(pastix)], 0 );
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

    isched_unsetaffinity();
    return NULL;
}

