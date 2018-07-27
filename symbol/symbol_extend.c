/**
 *
 * @file symbol_reorder.c
 *
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Vincent Bridonneau
 * @date 2018-07-16
 *
 * @precisions normal z -> s d c p
 *
 **/
#include "common.h"
#include "solver.h"
#include "symbol.h"
#include "symbol_reorder.h"
#include "pastix/order.h"
#include "blend/queue.h"
#include "blend/extendVector.h"

/* Find the column blok corresponding to the schur complement */
static inline symbol_cblk_t*
__find_cblk_schurfcol( symbol_matrix_t *symbptr )
{
    pastix_int_t   cblkid, cblknbr;
    symbol_cblk_t *cblk;

    cblk    = symbptr->cblktab;
    cblknbr = symbptr->cblknbr;
    for (cblkid = 0; cblkid <= cblknbr; ++cblkid)
    {
        if (cblk[cblkid].fcolnum == symbptr->schurfcol) {
            return cblk + cblkid;
        }
    }

    return NULL;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol
 *
 * @brief Extend the symbol matrix structure (compressed -> extented)
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The symbol structure to initialize.                                             *
 *******************************************************************************/
void
pastixSymbolExtend(       symbol_matrix_t *symbptr )
{
    /*Update block row and column indexes*/
    pastix_int_t   col, row, cblknbr, bloknbr, dof;
    symbol_cblk_t *cblk, *schurcblk;
    symbol_blok_t *blok;

    if ( symbptr == NULL ) {
        printf( "Warning : pastixSymbolExtend : symbmatrix is not initialised\n" );
        return;
    }

    pastixSymbolBase( symbptr, 0 );

    dof     = symbptr->dof;
    cblknbr = symbptr->cblknbr;
    cblk    = symbptr->cblktab;
    /* Update cblks first */
    for (col = 0; col <= cblknbr; ++col)
    {
        if (dof >= 1) {
            cblk[col].fcolnum *= dof;
            cblk[col].lcolnum *= dof;
        }
        else {
        }
    }

    bloknbr = symbptr->bloknbr;
    blok    = symbptr->bloktab;
    /* Update blocks then */
    for (row = 0; row < bloknbr; ++row)
    {
        if (dof >= 1) {
            blok[row].frownum *= dof;
            blok[row].lrownum *= dof;
        }
        else {
        }
    }

    if ( dof >= 1) {
        symbptr->nodenbr   *= dof;
        symbptr->schurfcol *= dof;
    }
    else {
        schurcblk = __find_cblk_schurfcol( symbptr );
    }

    symbptr->dof = 1;
}
