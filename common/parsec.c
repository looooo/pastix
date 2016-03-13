/**
 *
 * @file parsec.c
 *
 *  PaStiX PaRSECroutines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
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
#include <dague.h>
#include "dague/utils/mca_param.h"

void pastix_parsec_init( pastix_data_t *pastix,
                         int *argc, char **argv[] )
{
    pastix_int_t *iparm = pastix->iparm;
    char *value;
    int rc;

    if (iparm[IPARM_GPU_NBR] > 0) {
        rc = asprintf(&value, "%d", (int)(iparm[IPARM_GPU_NBR]));
        dague_setenv_mca_param( "device_cuda_enabled", value, &environ );

        if (iparm[IPARM_VERBOSE] > 2) {
            dague_setenv_mca_param( "device_show_statistics", "1", &environ );
        }
        if (iparm[IPARM_VERBOSE] > 3) {
            dague_setenv_mca_param( "device_show_capabilities", "1", &environ );
        }
    }
    pastix->parsec = dague_init( pastix->iparm[IPARM_THREAD_NBR],
                                 argc, argv );

    (void)rc;
}


void pastix_parsec_finalize( pastix_data_t *pastix )
{
    if (pastix->parsec != NULL) {
        dague_fini( &(pastix->parsec) );
    }
}
