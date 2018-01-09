/**
 *
 * @file order_apply_level_order.c
 *
 * PaStiX order function that apply reverse level ordering to the elimination
 * tree.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include <string.h>
#include "common.h"
#include "elimintree.h"
#include "pastix/order.h"

/**
 *******************************************************************************
 *
 * @ingroup order_dev
 *
 * @brief This routine build the elimination tree associated to an ordering.
 *
 *******************************************************************************
 *
 * @param[in] order
 *          The ordering structure.
 *
 * @param[out] nbroots
 *          On exit, contains the number of roots in the elimination tree
 *          (forest if > 1)
 *
 * @param[inout] roots
 *          On entry, array of size order->nodenbr.
 *          On exit, the first nbroots elements contains the roots of
 *          elimination trees.
 *
 *******************************************************************************
 *
 * @return the elimination tree structure associated to the ordering.
 *
 *******************************************************************************/
static inline EliminTree *
pastixOrderBuildEtree( const pastix_order_t *order,
                       pastix_int_t         *nbroots,
                       pastix_int_t         *roots )
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
 * @ingroup pastix_order
 *
 * @brief This routine reorder the elimination tree nodes per level.
 *
 *******************************************************************************
 *
 * @param[in] order
 *          The ordering structure to reorder.
 *
 * @param[in] level_tasks2d
 *          Define the ways 2D tasks are decided. If < 0, autolvel will be made
 *          based on all blocks above the minimal width_tasks2d criteria. If 0,
 *          1D tasks will be used, and if > 0, only the first level_tasks2d lvel
 *          of the elimination tree will be considered as 2D tasks.
 *
 * @param[in] width_tasks2d
 *          Define the minimal width for the supernodes that are considered as 2D
 *          blocks if level_tasks2d < 0. Unused otherwise.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if the ordering structure is incorrect.
 *
 *******************************************************************************/
int
pastixOrderApplyLevelOrder( pastix_order_t *order,
                            pastix_int_t    level_tasks2d,
                            pastix_int_t    width_tasks2d )
{
    pastix_order_t  oldorder;
    EliminTree     *etree;
    pastix_int_t    baseval;                  /* Node base value            */
    pastix_int_t    i, s, nbroots, node, sonsnbr;
    pastix_int_t    nfcol, ofcol, size;

    /* Parameter checks */
    if ( order == NULL ) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid order pointer");
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( (order->permtab == NULL) && (order->vertnbr > 0) ) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid order->permtab pointer");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( order->rangtab == NULL ) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid order->rangtab pointer");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( (order->treetab == NULL) && (order->cblknbr > 0) ) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid order->treetab pointer");
        return PASTIX_ERR_BADPARAMETER;
    }

    if (order->cblknbr < 0) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid nunber of column blocks");
        return PASTIX_ERR_BADPARAMETER;
    }
    baseval = order->baseval;
    if (baseval < 0) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid vertex node base number");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Quick return */
    if (order->cblknbr == 0) {
        return PASTIX_SUCCESS;
    }

    assert(baseval == order->rangtab[0]);

    memcpy( &oldorder, order, sizeof(pastix_order_t) );
    pastixOrderAlloc( order,
                      oldorder.vertnbr,
                      oldorder.cblknbr );

    /*
     * Build the elimination tree from top to bottom, and store the roots in the
     * permtab array
     */
    etree = pastixOrderBuildEtree( &oldorder,
                                   &nbroots,
                                   order->permtab );

    /*
     * Build the sorted array per level
     * If autolevel is enabled for 2D, we need to sort the 2D cblks first and
     * then the 1D.
     */
    if ( level_tasks2d < 0 )
    {
        pastix_int_t  pos_2D;
        pastix_int_t  pos_non_2D;
        pastix_int_t  sons2D;
        pastix_int_t  tot_nb_2D = 0;

        pastix_int_t *sorted = order->permtab;
        pastix_int_t *is_2D;
        MALLOC_INTERN(is_2D, order->cblknbr, pastix_int_t);
        memset(is_2D, 0, order->cblknbr * sizeof(pastix_int_t));

        for(i=0; i<nbroots; i++) {
            is_2D[order->cblknbr-1-i] = 1;
        }

        /* First pass to choose which nodes are 2D from the top to bottom */
        for(node=order->cblknbr-1; node>-1; node--) {

            /* Skip all nodes already set to 0 */
            while(is_2D[node] == 0 && node > 0){
                node--;
            }

            sonsnbr = etree->nodetab[node].sonsnbr;
            /*
             * If sonsnbr != 2, it is not a nested dissection node, and we
             * should not reorder them
             */
            if (sonsnbr == 2) {
                for(s=0; s<sonsnbr; s++) {
                    pastix_int_t son = eTreeSonI(etree, node, s);
                    size = oldorder.rangtab[ son+1 ] - oldorder.rangtab[ son ];
                    if (size >= width_tasks2d){
                        is_2D[son] = 1;
                        tot_nb_2D++;
                    }
                    else
                        is_2D[son] = 0;
                }
            }
        }

        pos_2D     = nbroots;
        pos_non_2D = tot_nb_2D+nbroots;

        /* Second pass to sort nodes: firstly by type (1D/2D) and then by levels */
        for(i=0; i<order->cblknbr; i++) {
            pastix_int_t current_2D     = 0;
            pastix_int_t current_non_2D = 0;

            node    = sorted[i];
            sonsnbr = etree->nodetab[node].sonsnbr;
            sons2D  = 0;

            size = oldorder.rangtab[ node+1 ] - oldorder.rangtab[ node ];

            for(s=0; s<sonsnbr; s++) {
                pastix_int_t son = eTreeSonI(etree, node, s);
                if (is_2D[son] == 1) {
                    sons2D++;
                }
            }

            /*
             * We put the sons in reverse order to keep the original order
             * betwen the brothers. This matters for the Minimum Degree part of
             * the ordering algorithm.
             */
            for(s=0; s<sonsnbr; s++) {
                pastix_int_t son = eTreeSonI(etree, node, s);
                if (is_2D[son] == 1){
                    sorted[pos_2D + sons2D - current_2D - 1] = son;
                    current_2D++;
                }
                else{
                    sorted[pos_non_2D + sonsnbr-sons2D - current_non_2D - 1] = son;
                    current_non_2D++;
                }
                etree->nodetab[ son ].fathnum = order->cblknbr - i - 1;
            }
            pos_2D     += sons2D;
            pos_non_2D += sonsnbr-sons2D;
        }
        memFree_null(is_2D);
    }
    else
    {
        pastix_int_t  pos = nbroots;
        pastix_int_t *sorted = order->permtab;

        for(i=0; i<order->cblknbr; i++) {
            node = sorted[i];
            sonsnbr = etree->nodetab[node].sonsnbr;
            /*
             * We put the sons in reverse order to keep the original order
             * betwen the brothers. This matters for the Minimum Degree part of
             * the ordering algorithm.
             */
            for(s=0; s<sonsnbr; s++) {
                pastix_int_t son = eTreeSonI(etree, node, s);
                sorted[ pos + sonsnbr-1-s ] = son;
                etree->nodetab[ son ].fathnum = order->cblknbr - i - 1;
            }
            pos += sonsnbr;
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

    pastixOrderExit( &oldorder );
    eTreeExit( etree );

    return PASTIX_SUCCESS;
}
