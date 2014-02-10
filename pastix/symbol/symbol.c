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
** $Id: symbol.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol.c                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the general purpose     **/
/**                routines for the symbolic matrix.       **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     07 oct 1998     **/
/**                # Version 0.1  : from : 03 dec 1998     **/
/**                                 to     03 dec 1998     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL

#include "common.h"
#include "symbol.h"

/******************************************/
/*                                        */
/* The symbolic matrix handling routines. */
/*                                        */
/******************************************/

/*+ This routine initializes the given
*** symbolic block matrix structure.
*** It returns:
*** - 0  : in all cases.
+*/

int
symbolInit (
SymbolMatrix * const        symbptr)
{
  memset (symbptr, 0, sizeof (SymbolMatrix));
#ifdef STARPU_GET_TASK_CTX
  symbptr->starpu_subtree_nbr=1;
#endif
  return (0);
}

/*+ This routine frees the contents
*** of the given symbolic block matrix.
*** It returns:
*** - VOID  : in all cases.
+*/

void
symbolExit (
SymbolMatrix * const        symbptr)
{
  if (symbptr->cblktab != NULL)
    memFree_null (symbptr->cblktab);
  if (symbptr->bloktab != NULL)
    memFree_null (symbptr->bloktab);

#ifdef SYMBOL_DEBUG
  symbolInit (symbptr);
#endif /* SYMBOL_DEBUG */
}

/*+ This routine reallocates the arrays
*** of the given symbolic block matrix.
*** It returns:
*** - VOID  : in all cases.
+*/

void
symbolRealloc (
SymbolMatrix * const        symbptr)
{
  SymbolCblk *        cblktab = NULL;
  SymbolBlok *        bloktab = NULL;

  if ((cblktab = (SymbolCblk *) memAlloc ((symbptr->cblknbr + 1) * sizeof (SymbolCblk))) == NULL)
    return;                                       /* Cannot move smallest array */
  memcpy  (cblktab, symbptr->cblktab, (symbptr->cblknbr + 1) * sizeof (SymbolCblk));
  memFree (symbptr->cblktab);                     /* Move column block array */
  symbptr->cblktab = cblktab;

  if ((bloktab = (SymbolBlok *) memAlloc (symbptr->bloknbr * sizeof (SymbolBlok))) == NULL)
    return;                                       /* Cannot move array */
  memcpy  (bloktab, symbptr->bloktab, symbptr->bloknbr * sizeof (SymbolBlok));
  memFree (symbptr->bloktab);                     /* Move column block array */
  symbptr->bloktab = bloktab;
}

/** Get face block for task E2 **/
pastix_int_t
symbolGetFacingBloknum(const SymbolMatrix *symbptr,
                       pastix_int_t bloksrc,
                       pastix_int_t bloknum,
                       pastix_int_t startsearch,
                       int ricar)
{
    pastix_int_t i;

    if(startsearch < symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum].bloknum )
        startsearch = symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum].bloknum;

    assert(startsearch < symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum+1].bloknum);

    if(ricar == 0)
    {
        for(i=startsearch;i<symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum+1].bloknum;i++)
            if(symbptr->bloktab[i].lrownum >= symbptr->bloktab[bloknum].frownum)
                break;

        assert( (symbptr->bloktab[i].frownum <= symbptr->bloktab[bloknum].frownum) &&
                (symbptr->bloktab[i].lrownum >= symbptr->bloktab[bloknum].lrownum) );

        return i;
    }
    else
    {
        for(i=startsearch;i<symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum+1].bloknum;i++)
        {
            if( (symbptr->bloktab[bloknum].frownum >= symbptr->bloktab[i].frownum && symbptr->bloktab[bloknum].frownum <= symbptr->bloktab[i].lrownum) ||
                (symbptr->bloktab[bloknum].lrownum >= symbptr->bloktab[i].frownum && symbptr->bloktab[bloknum].lrownum <= symbptr->bloktab[i].lrownum) ||
                (symbptr->bloktab[bloknum].frownum <= symbptr->bloktab[i].frownum && symbptr->bloktab[bloknum].lrownum >= symbptr->bloktab[i].lrownum) )
                return i;  /** We found the first block that matches **/
            if(symbptr->bloktab[bloknum].lrownum < symbptr->bloktab[i].frownum)
            {
                return -1;
            }
        }
    }
    return -1;
}

