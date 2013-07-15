/*
 *  File: order_compute_scotch.c
 *
 *  Wrapper to compute the ordering with Scotch Library.
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
#include <scotch.h>
#include "scotch_strats.h"
#include "csc_utils.h"
#include "cscd_utils_intern.h"


int orderComputeScotch( pastix_data_t *pastix_data, const pastix_graph_t *csc )
{
    Order        *ordemesh = pastix_data->ordemesh;
    SCOTCH_Graph  grafmesh;
    SCOTCH_Strat  stratdat;
    char          strat[1024];
    pastix_int_t *colptr;
    pastix_int_t *rows;
    pastix_int_t *iparm = pastix_data->iparm;
    pastix_int_t  procnum;
    pastix_int_t n;
    pastix_int_t nnz;
    int ret;

    procnum   = pastix_data->procnum;

    /* Check integer compatibility */
    if (sizeof(pastix_int_t) != sizeof(SCOTCH_Num)) {
        errorPrint("Inconsistent integer type\n");
        return INTEGER_TYPE_ERR;
    }

    n      = csc->n;
    colptr = csc->colptr;
    rows   = csc->rows;
    nnz    = colptr[n] - 1;

    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphInit <\n");
    orderInit(ordemesh, n, n);
    SCOTCH_graphInit( &grafmesh );

    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphBuild <\n");
    if (SCOTCH_graphBuild(&grafmesh,      /* Graph to build     */
                          1,              /* baseval            */
                          n,              /* Number of vertices */
                          colptr,         /* Vertex array       */
                          NULL,
                          NULL,           /* Array of vertex weights (DOFs) */
                          NULL,
                          nnz,            /* Number of arcs     */
                          rows,     /* Edge array         */
                          NULL))
        {
            errorPrint("pastix : graphBuildGraph");
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphCheck <\n");
    if (SCOTCH_graphCheck(&grafmesh)) {
        errorPrint("pastix: graphCheck");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }
    SCOTCH_graphBase(&grafmesh, 0);

    /* The graph is build, let's compute the ordering */
    SCOTCH_stratInit(&stratdat);

    /*
     * Create Strategy string for Scotch
     */
    /* default ordering */
    if (iparm[IPARM_DEFAULT_ORDERING] == API_YES) {
        if (iparm[IPARM_INCOMPLETE] == API_NO) {
            if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
                pastix_print(procnum, 0, "%s", "Scotch direct strategy\n");
            sprintf(strat, SCOTCH_STRAT_DIRECT);
        }
        else {
            if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
                pastix_print(procnum, 0, "%s", "Scotch incomplete strategy\n");
            sprintf(strat, SCOTCH_STRAT_INCOMP);
        }
    }
    /* personal ordering */
    else {
        sprintf(strat, SCOTCH_STRAT_PERSO,
                (long)  iparm[IPARM_ORDERING_SWITCH_LEVEL],
                (long)  iparm[IPARM_ORDERING_CMIN],
                (long)  iparm[IPARM_ORDERING_CMAX],
                ((float)iparm[IPARM_ORDERING_FRAT])/100,
                (long)  iparm[IPARM_ORDERING_SWITCH_LEVEL],
                (long)  iparm[IPARM_ORDERING_CMIN],
                (long)  iparm[IPARM_ORDERING_CMAX],
                ((float)iparm[IPARM_ORDERING_FRAT])/100);
        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
            pastix_print(procnum, 0, "Scotch personal strategy |%s|\n", strat);
    }

    ret = SCOTCH_stratGraphOrder (&stratdat, strat);
    if (ret == 0) {
        /* Compute graph ordering */
        ret = SCOTCH_graphOrderList(&grafmesh,
                                    (SCOTCH_Num)   n,
                                    (SCOTCH_Num *) NULL,
                                    &stratdat,
                                    (SCOTCH_Num *) ordemesh->permtab,
                                    (SCOTCH_Num *) ordemesh->peritab,
                                    (SCOTCH_Num *)&ordemesh->cblknbr,
                                    (SCOTCH_Num *) ordemesh->rangtab,
                                    NULL);

    }

    SCOTCH_stratExit (&stratdat);
    if (ret != 0) {           /* If something failed in Scotch */
        orderExit (ordemesh);    /* Free ordering arrays          */
        return INTERNAL_ERR;
    }

#if defined(FORGET_PARTITION)
    ordemesh->cblknbr = 0;
    if (ordemesh->rangtab != NULL) memFree_null(ordemesh->rangtab);
#else
    /* Redimensionnement de rangtab a cblknbr */
    ordemesh->rangtab =
        (pastix_int_t *) memRealloc (ordemesh->rangtab,
                                     (ordemesh->cblknbr + 1)*sizeof (pastix_int_t));
#endif

    return NO_ERR;
}
