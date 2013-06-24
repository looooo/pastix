/* Copyright INRIA 2004
**
** This file is part of the Scotch distribution.
**
** The Scotch distribution is libre/free software; you can
** redistribute it and/or modify it under the terms of the
** GNU Lesser General Public License as published by the
** Free Software Foundation; either version 2.1 of the
** License, or (at your option) any later version.
**
** The Scotch distribution is distributed in the hope that
** it will be useful, but WITHOUT ANY WARRANTY; without even
** the implied warranty of MERCHANTABILITY or FITNESS FOR A
** PARTICULAR PURPOSE. See the GNU Lesser General Public
** License for more details.
**
** You should have received a copy of the GNU Lesser General
** Public License along with the Scotch distribution; if not,
** write to the Free Software Foundation, Inc.,
** 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
**
** $Id: order.c 2 2004-06-02 14:05:03Z ramet $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : order.c                                 **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module computes orderings.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 20 aug 1998     **/
/**                                 to     24 sep 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "common.h"
#include "order.h"

/***********************************/
/*                                 */
/* The ordering handling routines. */
/*                                 */
/***********************************/

/* This routine initializes the given
** ordering structure.
** It returns:
** - 0  : in all cases.
*/

int orderInit ( Order * const ordeptr,
                pastix_int_t cblknbr,
                pastix_int_t vertnbr)
{
    memset(ordeptr, 0, sizeof(Order));

    MALLOC_INTERN(ordeptr->permtab, vertnbr,   pastix_int_t);
    MALLOC_INTERN(ordeptr->peritab, vertnbr,   pastix_int_t);
    MALLOC_INTERN(ordeptr->rangtab, cblknbr+1, pastix_int_t);

    return PASTIX_SUCCESS;
}

/* This routine frees the contents
** of the given ordering.
** It returns:
** - VOID  : in all cases.
*/

void orderExit (Order * const ordeptr)
{
    if (ordeptr->rangtab != NULL)
        memFree_null (ordeptr->rangtab);
    if (ordeptr->permtab != NULL)
        memFree_null (ordeptr->permtab);
    if (ordeptr->peritab != NULL)
        memFree_null (ordeptr->peritab);

#ifdef ORDER_DEBUG
    memSet (ordeptr, ~0, sizeof (Order));
#endif /* ORDER_DEBUG */
}
