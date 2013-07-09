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


int orderComputeScotch(pastix_data_t *pastix_data)
{
    Order        *ordemesh = &(pastix_data->ordemesh);
    SCOTCH_Graph *grafmesh = &(ordemesh->grafmesh);
    SCOTCH_Strat  stratdat;
    char          strat[1024];
    pastix_int_t *colptr, *colptr_schur;
    pastix_int_t *rows, *rows_schur;
    pastix_int_t *perm_schur, *invp_schur;
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

    n      = pastix_data->n2;
    colptr = colptr_schur = pastix_data->col2;
    rows   = rows_schur   = pastix_data->row2;
    nnz    = colptr[n] - 1;

    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphInit <\n");
    orderInit(ordemesh, n, n);

    if (iparm[IPARM_SCHUR]         == API_YES ||
        iparm[IPARM_ISOLATE_ZEROS] == API_YES) {

        /* Allocate separate pointers for shur complement */
        MALLOC_INTERN(colptr_schur, n+1,   pastix_int_t);
        MALLOC_INTERN(rows_schur,   nnz,   pastix_int_t);
        MALLOC_INTERN(perm_schur,   n,     pastix_int_t);
        MALLOC_INTERN(invp_schur,   n,     pastix_int_t);

        /* Backup colptr and rows */
        memcpy(colptr_schur, colptr, (n+1)*sizeof(pastix_int_t));
        memcpy(rows_schur,   rows,   nnz  *sizeof(pastix_int_t));

        CSC_isolate(n,
                    colptr_schur,
                    rows_schur,
                    pastix_data->nschur,
                    pastix_data->listschur,
                    perm_schur,
                    invp_schur);

        memFree_null(invp_schur);

        n   = n - pastix_data->nschur;
        nnz = colptr_schur[n] - 1;
    }

    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphBuild <\n");
    if (SCOTCH_graphBuild(grafmesh,       /* Graph to build     */
                          1,              /* baseval            */
                          n,              /* Number of vertices */
                          colptr_schur,   /* Vertex array       */
                          NULL,
                          NULL,           /* Array of vertex weights (DOFs) */
                          NULL,
                          nnz,            /* Number of arcs     */
                          rows_schur,     /* Edge array         */
                          NULL))
        {
            errorPrint("pastix : graphBuildGraph");
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphCheck <\n");
    if (SCOTCH_graphCheck(grafmesh)) {
        errorPrint("pastix: graphCheck");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }
    SCOTCH_graphBase(grafmesh, 0);

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
        ret = SCOTCH_graphOrderList(grafmesh,
                                    (SCOTCH_Num)   n,
                                    (SCOTCH_Num *) NULL,
                                    &stratdat,
                                    (SCOTCH_Num *) ordemesh->permtab,
                                    (SCOTCH_Num *) ordemesh->peritab,
                                    (SCOTCH_Num *)&ordemesh->cblknbr,
                                    (SCOTCH_Num *) ordemesh->rangtab,
                                    NULL);

        if (iparm[IPARM_SCHUR] == API_YES ||
            iparm[IPARM_ISOLATE_ZEROS] == API_YES) {

            pastix_int_t *tmpperm = NULL;
            pastix_int_t iter;

            /* Ajouter la permutation du Schur au permtab/peritab
               Et un bloc au rangtab.
            */
            assert( colptr_schur != colptr );
            assert( rows_schur != rows );
            memFree_null(colptr_schur);
            memFree_null(rows_schur);

            /* Restor n and nnz to global value */
            n   = pastix_data->n2;
            nnz = colptr[n] - 1;

            ordemesh->rangtab[ordemesh->cblknbr+1] = n;
            ordemesh->cblknbr++;

            for(iter = n-pastix_data->nschur; iter < n; iter++)
                ordemesh->permtab[iter] = iter;

#if !defined(NDEBUG)
            for(iter = 0; iter < n; iter++) {
                ASSERT(ordemesh->permtab[iter] < n,  MOD_SOPALIN);
                ASSERT(ordemesh->permtab[iter] > -1, MOD_SOPALIN);
                ASSERT(perm_schur[iter] < n,  MOD_SOPALIN);
                ASSERT(perm_schur[iter] > -1, MOD_SOPALIN);
            }
#endif

            /* Build permtab by composition of perm_shur and permtab */
            MALLOC_INTERN(tmpperm, n, pastix_int_t);
            for(iter = 0; iter < n; iter++)
                tmpperm[iter] = ordemesh->permtab[perm_schur[iter]];

            memcpy(ordemesh->permtab, tmpperm, n*sizeof(pastix_int_t));
            memFree_null(tmpperm);
            memFree_null(perm_schur);

            /* Build peritab by computing inverse of permtab */
            for(iter = 0; iter < n; iter++)
                ordemesh->peritab[ordemesh->permtab[iter]] = iter;

#if !defined(NDEBUG)
            for(iter = 0; iter < n; iter++) {
                ASSERT(ordemesh->peritab[iter] < n,  MOD_SOPALIN);
                ASSERT(ordemesh->peritab[iter] > -1, MOD_SOPALIN);
            }
#endif

            /* Need to rebuild the graph for fax */
            if (ordemesh->malgrf == 1) {
                SCOTCH_graphExit(grafmesh);
                ordemesh->malgrf = 0;
            }

            if (SCOTCH_graphBuild(grafmesh,       /* Graph to build     */
                                  1,              /* baseval            */
                                  n,              /* Number of vertices */
                                  colptr,         /* Vertex array       */
                                  NULL,
                                  NULL,           /* Array of vertex weights (DOFs) */
                                  NULL,
                                  nnz,            /* Number of arcs     */
                                  rows,           /* Edge array         */
                                  NULL))
                {
                    errorPrint("pastix : graphBuildGraph");
                    EXIT(MOD_SOPALIN,INTERNAL_ERR);
                }
            SCOTCH_graphBase(grafmesh, 0);
        }
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
