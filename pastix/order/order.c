/**
 *
 * @file order.c
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains functions to init/clean the order structure.
 *
 * @version 5.1.0
 * @author Francois Pellegrini
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "order.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering
 *
 * orderInit - Initialize the order structure. The base value is set to 0 by
 * default.
 *
 *******************************************************************************
 *
 * @param[in,out] ordeptr
 *          The data structure to initialize.
 *
 * @param[in] vertnbr
 *          The number of nodes, this is the size of the internal permtab and
 *          peritab arrays.
 *
 * @param[in] cblknbr
 *          The number of supernodes. The internal rangtab array is of size
 *          cblknbr+1.
 *          If cblknbr == 0, rangtab is not allocated.
 *
 *******************************************************************************/
int orderInit ( Order * const ordeptr,
                pastix_int_t vertnbr,
                pastix_int_t cblknbr)
{
    memset(ordeptr, 0, sizeof(Order));

    ordeptr->vertnbr = vertnbr;
    ordeptr->cblknbr = cblknbr;

    if (vertnbr != 0) {
        MALLOC_INTERN(ordeptr->permtab, vertnbr, pastix_int_t);
        MALLOC_INTERN(ordeptr->peritab, vertnbr, pastix_int_t);
    }

    if (cblknbr != 0) {
        MALLOC_INTERN(ordeptr->rangtab, cblknbr+1, pastix_int_t);
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering
 *
 * orderExit - Free the arrays initialized in the order structure.
 *
 *******************************************************************************
 *
 * @param[in,out] ordeptr
 *          The data structure to clean.
 *
 *******************************************************************************/
void orderExit (Order * const ordeptr)
{
    if (ordeptr->rangtab != NULL)
        memFree_null (ordeptr->rangtab);
    if (ordeptr->permtab != NULL)
        memFree_null (ordeptr->permtab);
    if (ordeptr->peritab != NULL)
        memFree_null (ordeptr->peritab);

    memset(ordeptr, 0, sizeof(Order) );
}
