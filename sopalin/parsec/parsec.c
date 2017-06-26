/**
 *
 * @file parsec.c
 *
 *  PaStiX PaRSEC routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 6.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @addtogroup pastix_parsec
 * @{
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "common.h"
#if !defined(PASTIX_WITH_PARSEC)
#error "This file should not be compiled if PaRSEC is not enabled"
#endif
#include <stdio.h>
#include <parsec.h>
#include "parsec/utils/mca_param.h"

/**
 *******************************************************************************
 *
 * @brief Startup the PaRSEC runtime system.
 *
 * This function initialize and startup the PaRSEC runtime system with PaStix
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
 ******************************************************************************/
void
pastix_parsec_init( pastix_data_t *pastix,
                    int *argc, char **argv[] )
{
    extern char **environ;
    pastix_int_t *iparm = pastix->iparm;
    char *value;
    int rc;

    if (iparm[IPARM_GPU_NBR] > 0) {
        rc = asprintf(&value, "%d", (int)(iparm[IPARM_GPU_NBR]));
        parsec_setenv_mca_param( "device_cuda_enabled", value, &environ );

        rc = asprintf(&value, "%d", (int)(iparm[IPARM_GPU_MEMORY_BLOCK_SIZE]));
        parsec_setenv_mca_param( "device_cuda_memory_block_size", value, &environ );

        rc = asprintf(&value, "%d", (int)(iparm[IPARM_GPU_MEMORY_PERCENTAGE]));
        parsec_setenv_mca_param( "device_cuda_memory_use", value, &environ );

        if (iparm[IPARM_VERBOSE] > 2) {
            parsec_setenv_mca_param( "device_show_statistics", "1", &environ );
        }
        if (iparm[IPARM_VERBOSE] > 3) {
            parsec_setenv_mca_param( "device_show_capabilities", "1", &environ );
        }
    }
    pastix->parsec = parsec_init( iparm[IPARM_THREAD_NBR],
                                  argc, argv );

    (void)rc;
}

/**
 *******************************************************************************
 *
 * @brief Finalize the PaRSEC runtime system.
 *
 * This function stop the PaRSEC runtime system.
 *
 *******************************************************************************
 *
 * @param[inout] pastix
 *          The main pastix_data structure.
 *
 ******************************************************************************/
void
pastix_parsec_finalize( pastix_data_t *pastix )
{
    if (pastix->parsec != NULL) {
        parsec_fini( (parsec_context_t**)&(pastix->parsec) );
    }
}

/**
 * @}
 */
