/**
 *
 * @file starpu.c
 *
 * @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2023-07-21
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#include "common.h"
#if !defined(PASTIX_WITH_STARPU)
#error "This file should not be compiled if Starpu is not enabled"
#endif
#include <stdio.h>
#include "pastix_starpu.h"

#if defined(PASTIX_STARPU_HETEROPRIO)
#include <starpu_heteroprio.h>
/**
 *******************************************************************************
 *
 * @brief Inits the heteroprio priorities, mappings and accelerations.
 *
 * This function initializes the heteroprio system in an arbitrary manner for now.
 * Used as a callback by starpu directly.
 *
 ******************************************************************************/
void
init_heteroprio( unsigned ctx )
{
    unsigned idx;
    /* CPU uses 4 buckets and visits them in the natural order */
    starpu_heteroprio_set_nb_prios( ctx, STARPU_CPU_WORKER, BucketNumber );
    /* It uses direct mapping idx => idx */
    for ( idx = 0; idx < BucketNumber; ++idx ) {
        starpu_heteroprio_set_mapping( ctx, STARPU_CPU_WORKER, idx, idx );
        /* If there are no CUDA workers, we must tell that CPU is faster */
        starpu_heteroprio_set_faster_arch( ctx, STARPU_CPU_WORKER, idx );
    }
    if ( starpu_cuda_worker_get_count() ) {
        const int   cuda_matching[] = { BucketGEMM2D, BucketTRSM2D, BucketGEMM1D };
        const float cuda_factor     = 125.0f / starpu_cpu_worker_get_count();
        const float cuda_factors[]  = { cuda_factor, cuda_factor, cuda_factor };
        /* CUDA is enabled and uses 2 buckets */
        starpu_heteroprio_set_nb_prios( ctx, STARPU_CUDA_WORKER, 3 );

        for ( idx = 0; idx < 3; ++idx ) {
            /* CUDA has its own mapping */
            starpu_heteroprio_set_mapping( ctx, STARPU_CUDA_WORKER, idx, cuda_matching[idx] );
            /* For its buckets, CUDA is the fastest */
            starpu_heteroprio_set_faster_arch( ctx, STARPU_CUDA_WORKER, cuda_matching[idx] );
            /* And CPU is slower by a factor dependant on the bucket */
            starpu_heteroprio_set_arch_slow_factor(
                ctx, STARPU_CPU_WORKER, cuda_matching[idx], cuda_factors[idx] );
        }
    }
}
#endif

/**
 *******************************************************************************
 *
 * @brief Startup the StarPU runtime system.
 *
 * This function initialize and startup the StarPU runtime system with PaStix
 * configuration variables
 *
 *******************************************************************************
 *
 * @param[inout] pastix
 *          The main pastix_data structure.
 *
 * @param[inout] argc
 *          The number of arguments of the main program.
 *
 * @param[inout] argv
 *          The list of argument given to the main program.
 *
 * @param[in] bindtab
 *          The binding array of size the number of threads if a specific
 *          binding is required, NULL otherwise.
 *
 ******************************************************************************/
void
pastix_starpu_init( pastix_data_t *pastix,
                    int *argc, char **argv[],
                    const int *bindtab )
{
    struct starpu_conf *conf;
    pastix_int_t *iparm = pastix->iparm;
    int rc;

    if ( pastix->starpu != NULL )
        return;

    pastix->starpu = malloc(sizeof(struct starpu_conf));
    starpu_conf_init( pastix->starpu );

#if defined(PASTIX_STARPU_PROFILING_LOG)
    pastix_gendirectories( pastix );
    profiling_log_init( pastix->dir_local );
#endif

    /* Force no GPUs if CUDA has not been enabled in PaStiX */
#if !defined(PASTIX_WITH_CUDA)
    iparm[IPARM_GPU_NBR] = 0;
#endif

    conf = pastix->starpu;
    conf->ncpus = pastix_imax( 1, (iparm[IPARM_THREAD_NBR] - iparm[IPARM_GPU_NBR] - 1) );
    conf->ncuda = iparm[IPARM_GPU_NBR];
    conf->nopencl = 0;

#if defined(PASTIX_STARPU_HETEROPRIO)
    /*
     * Set scheduling to heteroprio in any case if requested at compilation
     */
    conf->sched_policy_name = "heteroprio";
#if !defined(HAVE_STARPU_SCHED_POLICY_CALLBACK )
    conf->sched_policy_init = &init_heteroprio;
#else
    conf->sched_policy_callback = &init_heteroprio;
#endif
#else /* PASTIX_STARPU_HETEROPRIO */

    if (conf->ncuda > 0) {
#if defined(PASTIX_GENERATE_MODEL)
        pastix_print( pastix->procnum, 0,
                      "WARNING: PaStiX compiled with -DPASTIX_GENERATE_MODEL forces:\n"
                      "    - a single event per stream\n"
                      "    - a single stream per GPU\n"
                      "    - restore the automatic detection of the number of threads\n" );

        conf->ncpus = -1;

        pastix_setenv( "STARPU_NWORKER_PER_CUDA", "1", 1 );
        pastix_setenv( "STARPU_CUDA_PIPELINE", "0", 1 );
#endif

        conf->sched_policy_name = "dmdas";
    }
    else {
        /*
         * Set scheduling to "ws"/"lws" if no cuda devices used because it
         * behaves better on homogneneous architectures. If the user wants
         * to use another scheduling strategy, he can set STARPU_SCHED
         * env. var. to whatever he wants
         */
#if (STARPU_MAJOR_VERSION > 1) || ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION >= 2))
        conf->sched_policy_name = "lws";
#else
        conf->sched_policy_name = "ws";
#endif
    }
#endif /* PASTIX_STARPU_HETEROPRIO */

    if ( bindtab != NULL ) {
        int i;

        assert( iparm[IPARM_THREAD_NBR] < STARPU_NMAXWORKERS );
        conf->use_explicit_workers_bindid = 1;

        for(i=0; i < pastix_imin( iparm[IPARM_THREAD_NBR], STARPU_NMAXWORKERS ); i++) {
            conf->workers_bindid[i] = bindtab[i];
        }
    }
#if defined(STARPU_USE_FXT)
    starpu_fxt_autostart_profiling( 0 ); /* FxT starts profiling upon explicit call only */
#endif
    rc = starpu_init( conf );

    starpu_malloc_on_node_set_default_flags( STARPU_MAIN_RAM,
                                             STARPU_MALLOC_PINNED
                                             | STARPU_MALLOC_COUNT
#if defined(PASTIX_STARPU_SIMULATION)
                                             | STARPU_MALLOC_SIMULATION_FOLDED
#endif
                                             );

#if defined(PASTIX_WITH_MPI)
#if defined(PASTIX_DEBUG_MPI) && !defined(PASTIX_STARPU_SIMULATION)
    {
        int flag = 0;
        MPI_Initialized( &flag );
        assert( flag );
    }
#endif
    starpu_mpi_init_comm( argc, argv, 0, pastix->inter_node_comm );
    pastix_starpu_tag_init( pastix );
#endif
    pastix_starpu_interface_init();

#if defined(PASTIX_WITH_CUDA) && !defined(PASTIX_STARPU_SIMULATION)
    starpu_cublas_init();
#endif

    /* Suspend threads until we need them */
    starpu_pause();

    assert( pastix->starpu != NULL );

    (void)argc; (void)argv;
    (void)rc;
}

/**
 *******************************************************************************
 *
 * @brief Finalize the StarPU runtime system.
 *
 * This function stop the StarPU runtime system.
 *
 *******************************************************************************
 *
 * @param[inout] pastix
 *          The main pastix_data structure.
 *
 ******************************************************************************/
void
pastix_starpu_finalize( pastix_data_t *pastix )
{
    if ( pastix->starpu != NULL ) {
        starpu_resume();

        pastix_starpu_interface_fini();

        profiling_display_allinfo();
        profiling_log_fini();

#if defined(PASTIX_WITH_MPI)
        starpu_mpi_shutdown();
#endif
#if defined(PASTIX_WITH_CUDA) && !defined(PASTIX_STARPU_SIMULATION)
        starpu_cublas_shutdown();
#endif
        starpu_shutdown();

        free( pastix->starpu );
        pastix->starpu = NULL;
    }
}

/**
 * @}
 */
