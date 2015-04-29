/**
 *
 * @file symbol_check.c
 *
 *  Copyright Inria 1999-2015
 *
 *  PaStiX symbol structure routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Francois Pelegrin
 * @author David Goudin
 * @author Pascal Henon
 * @author Francois Pellegrini
 * @author Pierre Ramet
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
 * symbolCheck - Checks the consistency of the given symbolic block matrix.
 * Because of incomplete factorization, from version 1.0, no check is performed
 * regarding the existence of facing blocks in facing columns.
 *
 * TODO: Complete test set to check the brow informations.
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The symbol structure to check.
 *
 *******************************************************************************
 *
 * @return
 *          \retval 0 if the symbol matrix is correct
 *          \retval 1 if incorrect
 *
 *******************************************************************************/
int
symbolCheck(const SymbolMatrix * const  symbptr)
{
    pastix_int_t      baseval; /* Base value                           */
    const SymbolCblk *cblktax; /* Based access to cblktab              */
    pastix_int_t      cblkmax; /* Maximum column block index           */
    pastix_int_t      cblknum; /* Based number of current column block */
    const SymbolBlok *bloktax; /* Based access to bloktab              */
    pastix_int_t      blokmax; /* Maximum block index                  */
    pastix_int_t      bloknum; /* Based number of current block        */
    pastix_int_t      nodemax; /* Maximum node index                   */

    baseval = symbptr->baseval;
    cblktax = symbptr->cblktab - baseval;
    cblkmax = symbptr->cblknbr + (baseval - 1);
    bloktax = symbptr->bloktab - baseval;
    blokmax = symbptr->bloknbr + baseval;
    nodemax = symbptr->nodenbr;

    for (cblknum = bloknum = baseval;
         cblknum <= cblkmax; cblknum ++) {
        if ((cblktax[cblknum].fcolnum     <  baseval)                  ||
            (cblktax[cblknum].lcolnum     >  nodemax)                  ||
            (cblktax[cblknum].bloknum     >  blokmax)                  ||
            (cblktax[cblknum].fcolnum     >  cblktax[cblknum].lcolnum) ||
            (cblktax[cblknum + 1].brownum >= cblktax[cblknum].brownum) ||
            (cblktax[cblknum + 1].fcolnum <= cblktax[cblknum].lcolnum) ||
            (cblktax[cblknum + 1].bloknum <= cblktax[cblknum].bloknum)) {
            errorPrint ("symbolCheck: invalid column block array");
            assert(0);
            return     (1);
        }

        if ((bloktax[bloknum].frownum != cblktax[cblknum].fcolnum) ||
            (bloktax[bloknum].lrownum != cblktax[cblknum].lcolnum) ||
            (bloktax[bloknum].fcblknm != cblknum)) {
            errorPrint ("symbolCheck: invalid diagonal block");
            assert(0);
            return     (1);
        }

        for (bloknum ++; bloknum < cblktax[cblknum + 1].bloknum; bloknum ++) {
            if ((bloktax[bloknum].lcblknm != cblknum)                      ||
                (bloktax[bloknum].fcblknm <  baseval)                      ||
                (bloktax[bloknum].fcblknm >  cblkmax)                      ||
                (bloktax[bloknum].frownum <= bloktax[bloknum - 1].lrownum) ||
                (bloktax[bloknum].fcblknm <  bloktax[bloknum - 1].fcblknm)) {
                errorPrint ("symbolCheck: invalid block array");
                assert(0);
                return     (1);
            }
        }
    }

    assert( cblktax[cblknum].brownum == (symbptr->bloknbr - symbptr->cblknbr) );
    return (0);
}

