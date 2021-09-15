/**
 *
 * @file order_compute_personal.c
 *
 * PaStiX order driver to perform ordering with personal algorithm.
 *
 * @copyright 2004-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2021-09-15
 *
 **/
#include "common.h"
#include "graph/graph.h"
#include "pastix/order.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Computes the personal ordering of the PaStiX instance.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the field ordemesh is initialize with the result of the
 *          ordering.
 *
 * @param[inout] myorder
 *          On entry, the permutation provide by the user.
 *          On exit, if the structure attributs != NULL and IPARM_ORDERING parameter
 *          is not set to PastixOrderPersonal, contains the permutation generated.
 *          Otherwise, it is not referenced.
 *
 * @param[in] n
 *          The size of the ordemesh vertices.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_FILE if a problem occurs during the read in pastixOrderLoad
 *
 *******************************************************************************/
int
pastixOrderComputePersonal( pastix_data_t  *pastix_data,
                            pastix_order_t *myorder,
                            pastix_int_t    n )
{
    pastix_order_t *ordemesh;
    pastix_int_t   *iparm;
    pastix_int_t    i;
    int             procnum;
    int             retval = PASTIX_SUCCESS;

    ordemesh = pastix_data->ordemesh;
    iparm    = pastix_data->iparm;
    procnum  = pastix_data->procnum;

    /* Load ordering from file */
    if ( iparm[IPARM_IO_STRATEGY] & PastixIOLoad ) {
        if (iparm[IPARM_VERBOSE] > PastixVerboseNot) {
            pastix_print(procnum, 0, OUT_ORDER_METHOD, "Load" );
        }
        retval = pastixOrderLoad( pastix_data, ordemesh );
        return retval;
    }
    pastixOrderAlloc(ordemesh, n, 0);

    /* Rebase the Personal ordering to 0 */
    if ( myorder != NULL ) {
        assert( myorder != NULL );
        assert( myorder->vertnbr == n );
        pastixOrderBase(myorder, 0);
    }

    if ( (myorder == NULL) || (myorder->permtab == NULL) ) {
        if ( (myorder == NULL) || (myorder->peritab == NULL) ) {
            if (iparm[IPARM_VERBOSE] > PastixVerboseNot) {
                pastix_print(procnum, 0, OUT_ORDER_METHOD, "Personal (identity)" );
            }
            for(i=0; i<n; i++) {
                ordemesh->permtab[i] = i;
                ordemesh->peritab[i] = i;
            }
        }
        else {
            if (iparm[IPARM_VERBOSE] > PastixVerboseNot) {
                pastix_print(procnum, 0, OUT_ORDER_METHOD, "Personal (from myorder->peritab)" );
            }
            /* generate permtab from myorder->peritab */
            for(i=0;i<n;i++) {
                ordemesh->permtab[myorder->peritab[i]] = i;
            }
            memcpy(ordemesh->peritab, myorder->peritab, n*sizeof(pastix_int_t));
        }
    }
    else {
        if (myorder->peritab == NULL) {
            if (iparm[IPARM_VERBOSE] > PastixVerboseNot) {
                pastix_print(procnum, 0, OUT_ORDER_METHOD, "Personal (from myorder->permtab)" );
            }
            /* generate peritab from myorder->permtab */
            for(i=0;i<n;i++) {
                ordemesh->peritab[myorder->permtab[i]] = i;
            }
            memcpy(ordemesh->permtab, myorder->permtab, n*sizeof(pastix_int_t));
        }
        else {
            if (iparm[IPARM_VERBOSE] > PastixVerboseNot) {
                pastix_print(procnum, 0, OUT_ORDER_METHOD, "Personal (myorder->permtab/peritab)" );
            }
            memcpy(ordemesh->permtab, myorder->permtab, n*sizeof(pastix_int_t));
            memcpy(ordemesh->peritab, myorder->peritab, n*sizeof(pastix_int_t));
        }
    }

    /* Destroy the rangtab */
    ordemesh->cblknbr = 0;
    memFree_null( ordemesh->rangtab );
    /* Destroy the treetab */
    memFree_null( ordemesh->treetab );

    /* If treetab is provided, user must also provide rangtab */
    if ( myorder != NULL ) {
        assert( !( (myorder->rangtab == NULL) && (myorder->treetab != NULL) ) );
        if (myorder->rangtab != NULL )
        {
            ordemesh->cblknbr = myorder->cblknbr;
            MALLOC_INTERN(ordemesh->rangtab, myorder->cblknbr+1, pastix_int_t);
            memcpy(ordemesh->rangtab, myorder->rangtab, (myorder->cblknbr+1)*sizeof(pastix_int_t));
        }
        if (myorder->treetab != NULL )
        {
            MALLOC_INTERN(ordemesh->treetab, myorder->cblknbr, pastix_int_t);
            memcpy(ordemesh->treetab, myorder->treetab, myorder->cblknbr*sizeof(pastix_int_t));
        }
    }

    return retval;
}
