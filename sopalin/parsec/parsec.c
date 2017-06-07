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
 **/
#define _GNU_SOURCE
#include "common.h"
#if !defined(PASTIX_WITH_PARSEC)
#error "This file should not be compiled if PaRSEC is not enabled"
#endif
#include <stdio.h>
#include <parsec.h>
#include "parsec/utils/mca_param.h"

void pastix_parsec_init( pastix_data_t *pastix,
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


void pastix_parsec_finalize( pastix_data_t *pastix )
{
    if (pastix->parsec != NULL) {
        parsec_fini( (parsec_context_t**)&(pastix->parsec) );
    }
}
