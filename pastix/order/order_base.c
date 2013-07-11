/**
 *
 * @file order_base.c
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains function to adjust the base value of the ordering.
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
 * orderBase - This routine sets the base of the given ordering structure to the
 * given base value.
 *
 *******************************************************************************
 *
 * @param[in,out] ordeptr
 *          The ordering to rebase.
 *
 * @param[in] baseval
 *          The base value to be used.
 *
 *******************************************************************************/
void
orderBase (Order * const ordeptr,
	   pastix_int_t  baseval)
{
    pastix_int_t baseadj;                    /* Base adjust */
    pastix_int_t cblknum;
    pastix_int_t vertnum;

    baseadj = baseval - ordeptr->baseval; /* Set base adjust     */
    if (baseadj == 0)                     /* If base already set */
	return;

    if (ordeptr->rangtab != NULL) {
        for (cblknum = 0; cblknum <= ordeptr->cblknbr; cblknum ++)
            ordeptr->rangtab[cblknum] += baseadj;
    }

    if (ordeptr->permtab != NULL) {
	for (vertnum = 0; vertnum < ordeptr->vertnbr; vertnum ++)
	    ordeptr->permtab[vertnum] += baseadj;
    }
    if (ordeptr->peritab != NULL) {
	for (vertnum = 0; vertnum < ordeptr->vertnbr; vertnum ++)
	    ordeptr->peritab[vertnum] += baseadj;
    }
}
