/**
 *
 * @file bindthread.c
 *
 * Copyright (c) 2008-2014 The University of Bordeaux, IPB, LaBRI, Inria -
 *                         Bordeaux-Sud-Ouest.  All rights reserved.
 *
 * Copyright (c) 2010-2014 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *  PaStiX thread binding routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains basic functions to bind threads.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "bindthread.h"

#if defined(HAVE_HWLOC)
#include "pastix_hwloc.h"
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

int pastix_bindthread(int cpu, int ht)
{
#if defined(HAVE_HWLOC)
    {
        cpu = pastix_hwloc_bind_on_core_index(cpu, ht);
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
        tid_t self_ktid = thread_self ();
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


int pastix_bindthread_mask(hwloc_cpuset_t cpuset)
{
#if defined(HAVE_HWLOC)
    return pastix_hwloc_bind_on_mask_index(cpuset);
#else
    fprintf(stderr, "pastix_bindthread_mask withou HwLOC is not implemented yet\n");
#endif
}

