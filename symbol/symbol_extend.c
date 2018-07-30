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

/*
 * Find the column blok corresponding to the schur complement
 * return 1 if schurfcol == nodenbr, 0 otherwise
 */
static inline int
__find_cblk_schurfcol( symbol_cblk_t   *schurcblk,
                       symbol_matrix_t *symbptr )
{
    pastix_int_t   cblkid, cblknbr;
    symbol_cblk_t *cblk;

    if (symbptr->schurfcol == symbptr->nodenbr) {
        return 1;
    }

    cblk    = symbptr->cblktab;
    cblknbr = symbptr->cblknbr;
    for (cblkid = 0; cblkid <= cblknbr; ++cblkid)
    {
        if (cblk[cblkid].fcolnum == symbptr->schurfcol) {
            schurcblk = cblk + cblkid;
            return 0;
        }
    }

    schurcblk = NULL;
    assert ( symbptr->schurfcol == symbptr->nodenbr );

    return 0;
}

static inline void
__update_lcolnum( symbol_cblk_t *cblk,
                  pastix_int_t  *dofs,
                  pastix_int_t  *rangtab,
                  pastix_int_t  *peritab,
                  pastix_int_t   cblknum )
{
    pastix_int_t sum, i, begin, end;

    begin = rangtab[cblknum];
    end   = rangtab[cblknum + 1];
    sum   = 0;
    for (i = begin; i < end; ++i)
    {
        sum += dofs[peritab[i] + 1] - dofs[peritab[i]];
    }

    cblk->lcolnum = sum + cblk->fcolnum;
}

/* Should be call after update cblks */
static inline void
__update_bloks_in_cblk( symbol_blok_t *blok,     /* First blok in a cblk */
                        symbol_cblk_t *cblktabu, /* Updated */
                        symbol_cblk_t *cblktabp, /* Not updated */
                        pastix_int_t  *dofs,
                        pastix_int_t  *peritab )
{
    pastix_int_t   begin, end, sum, i;
    pastix_int_t   cblknum;
    symbol_cblk_t *fcblku, *fcblkp;
    pastix_int_t   rg_init, rg_end, j;

    cblknum = blok->lcblknm;
    begin   = cblktabp[cblknum].bloknum;
    end     = cblktabp[cblknum + 1].bloknum;

    for (i = begin; i < end; ++i, ++blok)
    {
        fcblku  = cblktabu + blok->fcblknm; /* Updated */
        fcblkp  = cblktabp + blok->fcblknm; /* Not updated */
        rg_init = fcblkp->fcolnum;
        rg_end  = blok->frownum;
        sum     = 0;

        /* First : compute the frownum */
        for (j = rg_init; j < rg_end; ++j)
        {
            sum += dofs[peritab[j] + 1] - dofs[peritab[j]];
        }
        blok->frownum = fcblku->fcolnum + sum;
        rg_end = blok->lrownum;

        /* Then : compute the lrownum */
        for (; j <= rg_end; ++j)
        {
            sum += dofs[peritab[j] + 1] - dofs[peritab[j]];
        }
        blok->lrownum = fcblku->fcolnum + sum;
    }
}

static inline void
__extend_dof_cte(       symbol_matrix_t *symbptr,
                  const pastix_order_t  *ordeptr )
{
    pastix_int_t   col, cblknbr, bloknbr, row, dof;
    symbol_cblk_t *cblk;
    symbol_blok_t *blok;

    dof     = symbptr->dof;
    cblk    = symbptr->cblktab;
    cblknbr = symbptr->cblknbr;

    /* Update cblks first */
    for (col = 0; col <= cblknbr; ++col)
    {
        cblk[col].fcolnum *= dof;
        cblk[col].lcolnum *= dof;
    }

    bloknbr = symbptr->bloknbr;
    blok    = symbptr->bloktab;
    /* Update block row and column indexes */
    for (row = 0; row < bloknbr; ++row)
    {
        blok[row].frownum *= dof;
        blok[row].lrownum *= dof;
    }

    symbptr->nodenbr   *= dof;
    symbptr->schurfcol *= dof;

    (void)ordeptr;
}

static inline void
__extend_dof_not_cte(       symbol_matrix_t *symbptr,
                      const pastix_order_t  *ordeptr )
{
    pastix_int_t   col, cblknbr;
    symbol_cblk_t *cblk, *schurcblk, *cblktab_cpy;
    symbol_blok_t *blok;
    pastix_int_t  *rangtab, *peritab, *dofs;
    pastix_int_t   is_schur_n;

    dofs    = symbptr->dofs;
    cblknbr = symbptr->cblknbr;
    cblk    = symbptr->cblktab;
    peritab = ordeptr->peritab;
    rangtab = ordeptr->rangtab;

    is_schur_n = __find_cblk_schurfcol( schurcblk, symbptr );

    /*
     * Make a copy of cblkatb to have both updated values (to compute (f/l)rows)
     * And previous ones (to have start indexes for bloks)
     */
    MALLOC_INTERN( cblktab_cpy, cblknbr + 1, symbol_cblk_t );
    memcpy ( cblktab_cpy, cblk, (cblknbr + 1) * sizeof(symbol_cblk_t) );

    /* Update cblks first */
    cblk[0].fcolnum = 0;
    __update_lcolnum( cblk, dofs, rangtab, peritab, 0 );

    for (col = 1; col < cblknbr; ++col)
    {
        cblk[col].fcolnum = cblk[col - 1].lcolnum + 1;
        __update_lcolnum( cblk + col, dofs, rangtab, peritab, col );
    }

    /* Update last column blok */
    cblk[col].fcolnum  = cblk[col - 1].lcolnum + 1;
    cblk[col].lcolnum  = cblk[col].fcolnum;

    blok    = symbptr->bloktab;
    /* Udpate bloks according to cblks */
    for (col = 0; col < cblknbr; ++col)
    {
        __update_bloks_in_cblk( blok + cblk[col].bloknum,
                                cblk, cblktab_cpy, dofs, peritab );
    }

    symbptr->nodenbr = cblk[symbptr->cblknbr].lcolnum;
    if ( is_schur_n ) {
        symbptr->schurfcol = symbptr->nodenbr;
    }
    else {
        symbptr->schurfcol = schurcblk->fcolnum;
    }

    memFree_null( cblktab_cpy );
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
 * @param[in] ordeptr
 *          The ordering structure providing rangtab for supernode
 *          and the permutation tab. This parameter must be given before
 *          its expension.
 *
 *******************************************************************************/
void
pastixSymbolExtend(       symbol_matrix_t *symbptr,
                    const pastix_order_t  *ordeptr )
{
    if ( symbptr == NULL ) {
        printf( "Warning : pastixSymbolExtend : symbmatrix is not initialised\n" );
        return;
    }

    pastixSymbolBase( symbptr, 0 );

    if ( symbptr->dof >= 1 ) {
        __extend_dof_cte( symbptr, ordeptr );
    }
    else {
        __extend_dof_not_cte( symbptr, ordeptr );
    }

    symbptr->dof = 1;
}
