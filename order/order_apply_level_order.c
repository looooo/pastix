/**
 *
 * @file order_apply_level_order.c
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains the function that apply reverse level ordering to the elimination
 * tree.
 *
 * @version 5.1.0
 * @author Francois Pellegrini
 * @date 2013-06-24
 *
 **/
#include <string.h>
#include "common.h"
#include "elimin.h"
#include "order.h"

static inline EliminTree *
orderBuildEtree( const Order *order,
                 pastix_int_t *nbroots,
                 pastix_int_t *roots )
{
    EliminTree *etree = NULL;
    eTreeNode_t *enode;
    pastix_int_t i, fathnum;
    pastix_int_t totalsonsnbr = 0;
    pastix_int_t sonstabcur;

    MALLOC_INTERN(etree, 1, EliminTree);
    eTreeInit(etree);

    etree->nodenbr = order->cblknbr;
    MALLOC_INTERN(etree->nodetab, etree->nodenbr, eTreeNode_t);
    memset( etree->nodetab, 0, etree->nodenbr * sizeof(eTreeNode_t) );
    enode = etree->nodetab;

    /* Initialize the structure fields */
    *nbroots = 0;
    for(i=0; i<order->cblknbr; i++, enode++)
    {
        fathnum = order->treetab[i];
        enode->fathnum = fathnum;
        enode->fsonnum = -1;
        if (fathnum != -1) {
            assert(fathnum < order->cblknbr );
            etree->nodetab[ fathnum ].sonsnbr ++;
            totalsonsnbr++;
        }
        else {
            roots[ *nbroots ] = i;
            (*nbroots)++;
        }
    }

    /* Check that we have at least one root */
    assert(totalsonsnbr == order->cblknbr-(*nbroots));

    if( totalsonsnbr > 0 ) {
        MALLOC_INTERN(etree->sonstab, totalsonsnbr, pastix_int_t);
    }

    /* Set the index of the first sons */
    sonstabcur = 0;
    for(i=0; i<order->cblknbr; i++)
    {
        etree->nodetab[i].fsonnum = sonstabcur;
        sonstabcur += etree->nodetab[i].sonsnbr;
    }
    assert(sonstabcur == totalsonsnbr);

    /* Fill the sonstab */
    /* No need to go to the root */
    for(i=0; i<order->cblknbr; i++)
    {
        fathnum = etree->nodetab[i].fathnum;
        if ( fathnum == -1 )
            continue;
        etree->sonstab[ etree->nodetab[ fathnum ].fsonnum++ ] = i;
    }

    /* Restore fsonnum fields */
    sonstabcur = 0;
    for(i=0; i<order->cblknbr; i++)
    {
        etree->nodetab[i].fsonnum = sonstabcur;
        sonstabcur += etree->nodetab[i].sonsnbr;
    }
    assert(sonstabcur == totalsonsnbr);

    return etree;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering
 *
 * orderApplyLevelOrder - This routine reorder the elimination tree nodes per
 * level.
 *
 *******************************************************************************
 *
 * @param[in] ordeptr
 *          The ordering structure to check.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCESS on successful exit.
 *          \retval PASTIX_ERR_BADPARAMETER if the ordering structure is incorrect.
 *
 *******************************************************************************/
int
orderApplyLevelOrder( Order *order )
{
    Order               oldorder;
    EliminTree         *etree;
    pastix_int_t        baseval;                  /* Node base value            */
    pastix_int_t        i, s, nbroots, node, sonsnbr;
    pastix_int_t        nfcol, ofcol, size;

    /* Parameter checks */
    if ( order == NULL ) {
        errorPrint ("orderCheck: invalid order pointer");
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( order->permtab == NULL ) {
        errorPrint ("orderCheck: invalid order->permtab pointer");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( order->rangtab == NULL ) {
        errorPrint ("orderCheck: invalid order->rangtab pointer");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( order->treetab == NULL ) {
        errorPrint ("orderCheck: invalid order->treetab pointer");
        return PASTIX_ERR_BADPARAMETER;
    }

    if (order->cblknbr < 0) {
        errorPrint ("orderCheck: invalid nunber of column blocks");
        return PASTIX_ERR_BADPARAMETER;
    }

    baseval = order->baseval;
    if (baseval < 0) {
        errorPrint ("orderCheck: invalid vertex node base number");
        return PASTIX_ERR_BADPARAMETER;
    }

    assert(baseval == order->rangtab[0]);

    memcpy( &oldorder, order, sizeof( Order ) );
    orderInit( order,
               oldorder.vertnbr,
               oldorder.cblknbr );

    /**
     * Build the elimination tree from top to bottom, and store the roots in the permatb array
     */
    etree = orderBuildEtree( &oldorder,
                             &nbroots,
                             order->permtab );

    /**
     * Build the sorted array per level
     */
    {
        pastix_int_t  pos = nbroots;
        pastix_int_t *sorted = order->permtab;

        for(i=0; i<order->cblknbr; i++) {
            node = sorted[i];
            sonsnbr = etree->nodetab[node].sonsnbr;
            for(s=0; s<sonsnbr; s++) {
                pastix_int_t son = eTreeSonI(etree, node, s);
                sorted[ pos ] = son;
                pos++;
                etree->nodetab[ son ].fathnum = order->cblknbr - i - 1;
            }
        }
        assert(pos == order->cblknbr);
    }

    /* Let's rebuild peritab, treetab, and rangtab */
    order->rangtab[0] = 0;

    for(i=0; i<order->cblknbr; i++ ) {
        node = order->permtab[order->cblknbr - i - 1];
        size = oldorder.rangtab[ node+1 ] - oldorder.rangtab[ node ];

        ofcol = oldorder.rangtab[ node ];
        nfcol = order->rangtab[ i ];
        order->rangtab[ i+1 ] = nfcol + size;
        order->treetab[ i ] = etree->nodetab[ node ].fathnum;

        memcpy( order->peritab + nfcol, oldorder.peritab + ofcol,
                size * sizeof( pastix_int_t ) );
    }

    /* Update the permutation */
    for (i=0; i<order->vertnbr; i++) {
        order->permtab[ order->peritab[i] ] = i;
    }

    orderExit( &oldorder );
    eTreeExit( etree );

    return PASTIX_SUCESS;
}
