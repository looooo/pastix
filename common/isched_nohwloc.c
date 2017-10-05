/**
 *
 * @file isched_nohwloc.c
 *
 * @copyright 2008-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2010-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 *
 * PaStiX thread binding routines without hwloc.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "isched.h"

#if defined(PASTIX_ARCH_COMPAQ)
#  include <sys/types.h>
#  include <sys/resource.h>
#  include <sys/processor.h>
#  include <sys/sysinfo.h>
#  include <machine/hal_sysinfo.h>
#  define X_INCLUDE_CXML
#elif defined(PASTIX_HAVE_SCHED_SETAFFINITY)
#  include <linux/unistd.h>
#  include <sched.h>
#elif defined(PASTIX_OS_MACOS)
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
#if defined(PASTIX_OS_LINUX) || defined(PASTIX_OS_FREEBSD) || defined(PASTIX_OS_AIX)

        sys_corenbr = sysconf(_SC_NPROCESSORS_ONLN);

#elif defined(PASTIX_OS_MACOS)

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
#elif defined(PASTIX_OS_WINDOWS)
        SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        sys_corenbr = sysinfo.dwNumberOfProcessors;
#endif
    }
    pthread_mutex_unlock(&mutextopo);

    return 0;
}

int isched_nohwloc_destroy(){
    return 0;
}

int isched_nohwloc_world_size()
{
    if ( !topo_initialized )
        isched_nohwloc_init();
    return sys_corenbr;
}

int isched_nohwloc_bind_on_core_index(int cpu)
{
    if( -1 == cpu )  /* Don't try binding if not required */
        return -1;

    if ( !topo_initialized )
        isched_nohwloc_init();

#if defined(PASTIX_HAVE_SCHED_SETAFFINITY)
    {
        cpu_set_t mask;
        CPU_ZERO(&mask);
        CPU_SET(cpu, &mask);

#if defined(PASTIX_HAVE_OLD_SCHED_SETAFFINITY)
        if(sched_setaffinity(0,&mask) < 0)
#else /* PASTIX_HAVE_OLD_SCHED_SETAFFINITY */
        if(sched_setaffinity(0,sizeof(mask),&mask) < 0)
#endif /* PASTIX_HAVE_OLD_SCHED_SETAFFINITY */
        {
            return -1;
        }
    }
#elif defined(PASTIX_ARCH_PPC)
    {
        tid_t self_ktid = thread_self();
        bindprocessor(BINDTHREAD, self_ktid, cpu*2);
    }
#elif defined(PASTIX_ARCH_COMPAQ)
    {
        bind_to_cpu_id(getpid(), cpu, 0);
    }
#elif defined(PASTIX_OS_MACOS)
    {
        thread_affinity_policy_data_t ap;
        int                           ret;

        ap.affinity_tag = 1; /* non-null affinity tag */
        ret = thread_policy_set(
            mach_thread_self(),
            THREAD_AFFINITY_POLICY,
            (integer_t*) &ap,
            THREAD_AFFINITY_POLICY_COUNT );
        if(ret != 0) {
            return -1;
        }
    }
#endif /* Architectures */

    return cpu;
}

int isched_nohwloc_unbind()
{
    if ( !topo_initialized )
        isched_nohwloc_init();

#if !defined(PASTIX_AFFINITY_DISABLE)
#if defined(PASTIX_OS_LINUX) || defined(PASTIX_OS_FREEBSD)
    {
        int i;
#if defined(PASTIX_OS_LINUX)
        cpu_set_t set;
#elif defined(PASTIX_OS_FREEBSD)
        cpuset_t set;
#endif
        CPU_ZERO( &set );

        for(i=0; i<sys_corenbr; i++)
            CPU_SET( i, &set );

#if defined(HAVE_OLD_SCHED_SETAFFINITY)
        if( sched_setaffinity( 0, &set ) < 0 )
#else /* HAVE_OLD_SCHED_SETAFFINITY */
#if defined(PASTIX_OS_LINUX)
        if( sched_setaffinity( 0, sizeof(set), &set) < 0 )
#elif defined(PASTIX_OS_FREEBSD)
        if( cpuset_setaffinity(CPU_LEVEL_WHICH, CPU_WHICH_PID, 0, sizeof(set), &set) < 0 )
#endif
#endif /* HAVE_OLD_SCHED_SETAFFINITY */
            {
                pastix_warning("pastix_unsetaffinity", "Could not unbind thread");
                return -1;
            }

        return 0;
    }
#elif defined(PASTIX_OS_MACOS)
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
        /*     pastix_warning("pastix_unsetaffinity", "Could not unbind thread"); */
        /*     return PASTIX_ERR_UNKNOWN; */
        /* } */

        return 0;
    }
#elif defined(PASTIX_OS_WINDOWS)
    {
        int i;
        DWORD mask = 0;

        for(i=0; i<sys_corenbr; i++)
            mask |= 1 << i;

        if( SetThreadAffinityMask(GetCurrentThread(), mask) == 0) {
            pastix_warning("pastix_unsetaffinity", "Could not unbind thread");
            return -1;
        }
        return 0;
    }
#elif defined(PASTIX_OS_AIX)
    {
        /* TODO: check how to unbind the main thread if necessary for OpenMP */
        /* tid_t self_ktid = thread_self (); */
        /* bindprocessor(BINDTHREAD, self_ktid, rank); */
        return 0;
    }
#else
    return -1;
#endif
#endif /* !defined(PASTIX_AFFINITY_DISABLE) */

    return 0;
}
