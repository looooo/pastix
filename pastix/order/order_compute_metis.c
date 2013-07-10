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
#include <metis.h>

int orderComputeMetis( pastix_data_t *pastix_data, pastix_csc_t *csc )
{
    pastix_int_t *iparm    =   pastix_data->iparm;
    Order        *ordemesh = &(pastix_data->ordemesh);
    pastix_int_t  itervert;
    idx_t baseval = csc->colptr[0];
    idx_t opt[8];

    if ( sizeof(pastix_int_t) != sizeof(idx_t)) {
        errorPrint("Inconsistent integer type between PaStiX and Metis\n");
        return INTEGER_TYPE_ERR;
    }

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        pastix_print(procnum, 0, "%s", "Ordering: calling metis...\n");

    /* call METIS and fill ordemesh (provide a partition) */
    opt[OPTION_PTYPE  ] = (iparm[IPARM_DEFAULT_ORDERING] == API_YES) ? 0 : 1;

    /* TODO: Check without this first line that reset to 0 if default */
    opt[OPTION_PTYPE  ] = 0;
    opt[OPTION_CTYPE  ] = iparm[IPARM_ORDERING_SWITCH_LEVEL];
    opt[OPTION_ITYPE  ] = iparm[IPARM_ORDERING_CMIN];
    opt[OPTION_RTYPE  ] = iparm[IPARM_ORDERING_CMAX];
    opt[OPTION_DBGLVL ] = iparm[IPARM_ORDERING_FRAT];
    opt[OPTION_OFLAGS ] = iparm[IPARM_STATIC_PIVOTING];
    opt[OPTION_PFACTOR] = iparm[IPARM_METIS_PFACTOR];
    opt[OPTION_NSEPS  ] = iparm[IPARM_NNZEROS];

    /*METIS_NodeND(&n,verttab,edgetab,&baseval,opt,
      ordemesh->permtab,ordemesh->peritab);*/
    METIS_NodeND( &n, csc->colptr, csc->rows, &baseval,
                  opt, ordemesh->peritab, ordemesh->permtab);

    for (itervert=0; itervert<n+1; itervert++)
        ordemesh->rangtab[itervert] = itervert;
    ordemesh->cblknbr = n;

    return NO_ERR;
}
#endif /* defined(HAVE_METIS) */
