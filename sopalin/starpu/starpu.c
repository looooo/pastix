/**
 *
 * @file starpu.c
 *
 *  PaStiX StarPU routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 6.0.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
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
#include <starpu.h>

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

    conf = pastix->starpu;
    conf->ncpus = iparm[IPARM_THREAD_NBR];
    conf->ncuda = iparm[IPARM_GPU_NBR];
    conf->nopencl = 0;

    if (conf->ncuda > 0) {
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

    if ( bindtab != NULL ) {
        int i;

        assert( iparm[IPARM_THREAD_NBR] < STARPU_NMAXWORKERS );
        conf->use_explicit_workers_bindid = 1;

        for(i=0; i < pastix_imin( iparm[IPARM_THREAD_NBR], STARPU_NMAXWORKERS ); i++) {
            conf->workers_bindid[i] = bindtab[i];
        }
    }
    rc = starpu_init( conf );

    starpu_malloc_on_node_set_default_flags( STARPU_MAIN_RAM,
                                             STARPU_MALLOC_PINNED
                                             | STARPU_MALLOC_COUNT
#if defined(PASTIX_STARPU_SIMULATION)
                                             | STARPU_MALLOC_SIMULATION_FOLDED
#endif
                                             );

#if defined(PASTIX_WITH_MPI)
    {
        int flag = 0;
#if !defined(PASTIX_STARPU_SIMULATION)
        MPI_Initialized( &flag );
#endif
        starpu_mpi_init( argc, argv, !flag );
    }
#endif

#if defined(PASTIX_WITH_CUDA) && !defined(PASTIX_STARPU_SIMULATION)
    starpu_cublas_init();
#endif

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
    if (pastix->starpu != NULL) {
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
