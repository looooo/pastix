/*
 *  File: order_compute_metis.c
 *
 *  Function that handle the ordering step of PaStiX. It is a wrapper
 *  to the different ordering library that can be used by PaStiX.
 *
 *  Authors:
 *    Mathieu  Faverge    - faverge@labri.fr
 *    Xavier   LACOSTE    - lacoste@labri.fr
 *    Pierre   RAMET      - ramet@labri.fr
 *
 *  Dates:
 *    Version 0.0 - from 08 may 1998
 *                  to   08 jan 2001
 *    Version 1.0 - from 06 jun 2002
 *                  to   06 jun 2002
 */
#include "common.h"
#include "graph.h"
#include "order.h"
#include <metis.h>

int orderComputeMetis( pastix_data_t *pastix_data, const pastix_graph_t *graph )
{
    pastix_int_t *iparm    = pastix_data->iparm;
    Order        *ordemesh = pastix_data->ordemesh;
    pastix_int_t  procnum  = pastix_data->procnum;
    pastix_int_t  n;
    pastix_int_t  baseval = graph->colptr[0];
    idx_t opt[METIS_NOPTIONS];
    int rc;

    if ( sizeof(pastix_int_t) != sizeof(idx_t)) {
        errorPrint("Inconsistent integer type between PaStiX and Metis\n");
        return PASTIX_ERR_INTEGER_TYPE;
    }

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        pastix_print(procnum, 0, "%s", "Ordering: calling metis...\n");

    /* Former options used with Metis < 5.0.0 */
    /* opt[METIS_OPTION_CTYPE  ] = iparm[IPARM_ORDERING_SWITCH_LEVEL]; */
    /* opt[METIS_OPTION_ITYPE  ] = iparm[IPARM_ORDERING_CMIN]; */
    /* opt[METIS_OPTION_RTYPE  ] = iparm[IPARM_ORDERING_CMAX]; */
    /* opt[METIS_OPTION_DBGLVL ] = iparm[IPARM_ORDERING_FRAT]; */
    /* opt[METIS_OPTION_OFLAGS ] = iparm[IPARM_STATIC_PIVOTING]; */
    /* opt[METIS_OPTION_PFACTOR] = iparm[IPARM_METIS_PFACTOR]; */
    /* opt[METIS_OPTION_NSEPS  ] = iparm[IPARM_NNZEROS]; */

    /* Set of valid options for METIS_NodeND */
    METIS_SetDefaultOptions(opt);
    if (iparm[IPARM_DEFAULT_ORDERING] != API_YES) {
        /* opt[METIS_OPTION_CTYPE    ] = METIS_CTYPE_RM; /\* METIS_CTYPE_RM or METIS_CTYPE_SHEM                *\/ */
        /* opt[METIS_OPTION_RTYPE    ] = METIS_RTYPE_FM; /\* METIS_RTYPE_FM, _GREEDY, _SEP2SIDED or _SEP1SIDED *\/ */
        /* opt[METIS_OPTION_NO2HOP   ] = 0; */
        /* opt[METIS_OPTION_NSEPS    ] = 1;  /\* Default: 1  *\/ */
        /* opt[METIS_OPTION_NITER    ] = 10; /\* Default: 10 *\/ */
        /* opt[METIS_OPTION_UFACTOR  ] = ;  */
        /* opt[METIS_OPTION_COMPRESS ] = 0; */
        /* opt[METIS_OPTION_CCORDER  ] = 1; */
        opt[METIS_OPTION_PFACTOR  ] = iparm[IPARM_METIS_PFACTOR];
    }
    opt[METIS_OPTION_SEED     ] = 3452;
    opt[METIS_OPTION_NUMBERING] = baseval;
    opt[METIS_OPTION_DBGLVL   ] = 0;

    n = graph->n;
    orderInit( ordemesh, graph->n, 0 );
    ordemesh->baseval = baseval;
    rc = METIS_NodeND( &n, graph->colptr, graph->rows, NULL,
                       opt, ordemesh->peritab, ordemesh->permtab);

    assert( n == graph->n );
    if (rc != METIS_OK )
    {
        errorPrint("orderComputeMetis: Invalid code returned by METIS_NodeND (%d)\n", rc);
        return PASTIX_ERR_BADPARAMETER;
    }

#if defined(PASTIX_DEBUG_ORDERING)
    {
        pastix_int_t i;
        for(i=0; i<n; i++) {
            assert( ordemesh->permtab[i] >= baseval );
            assert( ordemesh->permtab[i] <  (n+baseval) );
            assert( ordemesh->peritab[i] >= baseval );
            assert( ordemesh->peritab[i] <  (n+baseval) );
        }
    }
#endif
    return PASTIX_SUCCESS;
}
