/**
 *
 * @file symbol_base.c
 *
 *  Copyright Inria 1999-2015
 *
 *  PaStiX symbol structure routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Francois Pelegrin
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "symbol.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol
 *
 * symbolBase - Sets the base of the given symbol matrix structure to the given
 * base value.
 *
 *******************************************************************************
 *
 * @param[in,out] symbptr
 *          The symbol structure to initialize.
 *
 * @param[in] baseval
 *          The base value.
 *
 *******************************************************************************/
void
symbolBase ( SymbolMatrix * const symbptr,
             const pastix_int_t   baseval)
{
    SymbolCblk  *cblk;
    SymbolBlok  *blok;
    pastix_int_t baseadj; /* Base adjust */
    pastix_int_t cblknum;
    pastix_int_t bloknum;

    baseadj = baseval - symbptr->baseval; /* Set base adjust     */
    if (baseadj == 0)                     /* If base already set */
        return;

    symbptr->baseval = baseval;           /* Set graph base */

    cblk = symbptr->cblktab;
    for (cblknum = 0; cblknum <= symbptr->cblknbr; cblknum ++, cblk++) {
        cblk->fcolnum += baseadj;
        cblk->lcolnum += baseadj;
        cblk->bloknum += baseadj;
    }

    blok = symbptr->bloktab;
    for (bloknum = 0; bloknum < symbptr->bloknbr; bloknum ++) {
        blok->frownum += baseadj;
        blok->lrownum += baseadj;
        blok->lcblknm += baseadj;
        blok->fcblknm += baseadj;
    }
}
