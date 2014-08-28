/**
 *
 * @file pastix_task_sopalin.c
 *
 *  PaStiX analyse routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "csc.h"
#include "bcsc.h"

/*
  Function: pastix_task_sopalin

  Distribution task.

  Parameters:
  pastix_data - PaStiX data structure
  pastix_comm - PaStiX MPI communicator

*/
void pastix_task_sopalin( pastix_data_t *pastix_data,
                          pastix_csc_t  *csc )
{
/* #ifdef PASTIX_WITH_MPI */
/*     MPI_Comm       pastix_comm = pastix_data->inter_node_comm; */
/* #endif */
/*     pastix_int_t   procnum  = pastix_data->inter_node_procnum; */
    pastix_int_t  *iparm    = pastix_data->iparm;
/*     double        *dparm    = pastix_data->dparm; */
/*     SolverMatrix  *solvmatr = pastix_data->solvmatr; */

/*     if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO) */
/*     { */
/*         switch(iparm[IPARM_FACTORIZATION]) */
/*         { */
/*         case API_FACT_LU: */
/*             print_onempi("%s", OUT_STEP_NUMFACT_LU); */
/*             break; */
/*         case API_FACT_LLT: */
/*             print_onempi("%s", OUT_STEP_NUMFACT_LLT); */
/*             break; */
/*         case API_FACT_LDLH: */
/*             print_onempi("%s", OUT_STEP_NUMFACT_LDLH); */
/*             break; */
/*         case API_FACT_LDLT: */
/*         default: */
/*             print_onempi("%s", OUT_STEP_NUMFACT_LDLT); */
/*         } */
/*     } */

    /**
     * Fill in the internal blocked CSC. We consider that if this step is called
     * the csc values have changed so we need to update the blocked csc.
     */
    if (pastix_data->bcsc != NULL)
    {
        bcscExit( pastix_data->bcsc );
        memFree_null( pastix_data->bcsc );
    }

    MALLOC_INTERN( pastix_data->bcsc, 1, pastix_bcsc_t );

    /**
     * Check that IPARM_ONLY_RAFF is not set, because we don't need to
     * initialize the A^t array in this case and we should have never arrived
     * here.
     */
    assert( !iparm[IPARM_ONLY_RAFF] );

    bcscInit( csc,
              pastix_data->ordemesh,
              pastix_data->solvmatr,
              1,
              pastix_data->bcsc );

    if ( iparm[IPARM_FREE_CSCUSER] ) {
        //TODO: cscExit( csc );
    }

    iparm[IPARM_START_TASK]++;
}
