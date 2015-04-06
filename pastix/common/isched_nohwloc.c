/**
 *
 * @file isched_nohwloc.c
 *
 * Copyright (c) 2008-2015 The University of Bordeaux, IPB, LaBRI, Inria -
 *                         Bordeaux-Sud-Ouest.  All rights reserved.
 *
 * Copyright (c) 2010-2015 The University of Tennessee and The University
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

#if defined(ARCH_COMPAQ)
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
#endif  /* define(ARCH_COMPAQ) */


static pthread_mutex_t  mutextopo = PTHREAD_MUTEX_INITIALIZER;
static volatile int sys_corenbr = 1;
static volatile int topo_initialized = 0;

int isched_nohwloc_init() {
    pthread_mutex_lock(&mutextopo);
    if ( !topo_initialized ) {
#if (defined PLASMA_OS_LINUX) || (defined PLASMA_OS_FREEBSD) || (defined PLASMA_OS_AIX)

        sys_corenbr = sysconf(_SC_NPROCESSORS_ONLN);

#elif (defined PLASMA_OS_MACOS)

        int mib[4];
        int cpu;
        size_t len = sizeof(cpu);

        /* set the mib for hw.ncpu */
        mib[0] = CTL_HW;
        mib[1] = HW_AVAILCPU;

        /* get the number of CPUs from the system */
        sysctl(mib, 2, &cpu, &len, NULL, 0);
        if( cpu < 1 ) {
            mib[1] = HW_NCPU;
            sysctl( mib, 2, &cpu, &len, NULL, 0 );
        }
        if( cpu < 1 ) {
            cpu = 1;
        }
        sys_corenbr = cpu;
#elif (defined PLASMA_OS_WINDOWS)
        SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        sys_corenbr = sysinfo.dwNumberOfProcessors;
#endif
    }
    pthread_mutex_unlock(&mutextopo);
}

void isched_nohwloc_finalize(){
    return;
}

int isched_nohwloc_bind_on_core_index(int cpu)
{
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

    return cpu;
}

int isched_nohwloc_unbind()
{
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
}
