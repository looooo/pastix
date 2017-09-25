/**
 *
 * @file elimintree.c
 *
 *  PaStiX analyse routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains basic functions to manipulate elimination tree structure.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @addtogroup blend_dev_elim
 * @{
 *
 **/
#include "common.h"
#include "symbol.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"

static inline pastix_int_t
compress_getNodeNbr( const EliminTree *etree,
                            const Cand       *candtab,
                            pastix_int_t      rootnum )
{
    pastix_int_t i, merge, fcand, lcand;
    pastix_int_t sonsnbr, nbnodes;

    fcand = candtab[rootnum].fcandnum;
    lcand = candtab[rootnum].lcandnum;

    sonsnbr = etree->nodetab[rootnum].sonsnbr;
    merge   = 1;
    nbnodes = 1;
    for( i=0; i<sonsnbr; i++ )
    {
        pastix_int_t son = eTreeSonI(etree, rootnum, i);
        nbnodes += compress_getNodeNbr( etree, candtab, son );

        if ( (fcand != candtab[son].fcandnum) ||
             (lcand != candtab[son].lcandnum) )
        {
            merge = 0;
        }
    }

    /* If all sons have the same candidate set, then they will be merged within the current node */
    if (merge) {
        nbnodes -= sonsnbr;
    }
    return nbnodes;
}

static inline void
compress_setSonsNbr( const EliminTree *etree,
                     const Cand       *candtab,
                     pastix_int_t      rootnum,
                     EliminTree       *ctree,
                     Cand             *ccand,
                     pastix_int_t      fathnum,
                     pastix_int_t     *cnodeidx,
                     pastix_int_t     *tmp )
{
    eTreeNode_t *cnode;
    double total;
    pastix_int_t i, merge, fcand, lcand;
    pastix_int_t sonsnbr, gdsonsnbr;
    (*cnodeidx)++;
    assert( *cnodeidx < ctree->nodenbr );

    fcand = candtab[rootnum].fcandnum;
    lcand = candtab[rootnum].lcandnum;

    /* Copy current node */
    cnode = ctree->nodetab + (*cnodeidx);
    cnode->total   = etree->nodetab[rootnum].total;
    cnode->subtree = etree->nodetab[rootnum].subtree;
    cnode->cripath = etree->nodetab[rootnum].cripath;
    cnode->fathnum = fathnum;

    ccand[ *cnodeidx ].fcandnum = fcand;
    ccand[ *cnodeidx ].lcandnum = lcand;

    sonsnbr = etree->nodetab[rootnum].sonsnbr;
    if (sonsnbr == 0) {
        return;
    }

    memcpy( tmp, etree->sonstab + etree->nodetab[rootnum].fsonnum,
            sonsnbr * sizeof(pastix_int_t) );
    do {
        total     = 0.;
        gdsonsnbr = 0;
        merge     = 1;
        for( i=0; i<sonsnbr; i++ )
        {
            pastix_int_t son = tmp[i];

            /* Backup grandsons after the sons */
            memcpy( tmp + sonsnbr + gdsonsnbr,
                    etree->sonstab + etree->nodetab[ son ].fsonnum,
                    etree->nodetab[ son ].sonsnbr * sizeof(pastix_int_t) );
            gdsonsnbr += etree->nodetab[son].sonsnbr;
            total     += etree->nodetab[son].total;

            if ( (fcand != candtab[son].fcandnum) ||
                 (lcand != candtab[son].lcandnum) )
            {
                merge = 0;
            }
        }

        /* If all sons have the same candidate set, then they will be merged within the current node */
        if ( merge ) {
            /* Grandsons become sons */
            for (i=0; i<gdsonsnbr; i++) {
                tmp[i] = tmp[sonsnbr+i];

            }
            sonsnbr = gdsonsnbr;
            cnode->total += total;
        }
    }
    while( merge && (sonsnbr>0) );

    /* Recurse on sons */
    merge = *cnodeidx;
    for(i=0; i<sonsnbr; i++) {
        pastix_int_t son = (*cnodeidx) + 1;
        compress_setSonsNbr( etree, candtab, tmp[i],
                             ctree, ccand,
                             merge, cnodeidx, tmp+sonsnbr );
        tmp[i] = son;
    }
    cnode->sonsnbr = sonsnbr;

    /* Compress the single candidate nodes */
    for(i=0; i<sonsnbr; i++) {
        pastix_int_t j, soni, sonj;

        soni = tmp[i];
        if ( ccand[soni].fcandnum != ccand[soni].lcandnum )
            continue;

        for( j=i+1; j<sonsnbr; j++ ) {
            sonj = tmp[j];

            if ( (ccand[sonj].fcandnum != ccand[soni].fcandnum) ||
                 (ccand[sonj].lcandnum != ccand[soni].lcandnum) )
                continue;

            ctree->nodetab[soni].total   += ctree->nodetab[sonj].total;
            ctree->nodetab[soni].subtree += ctree->nodetab[sonj].subtree;
            assert( ctree->nodetab[sonj].sonsnbr == 0 );
            assert( ctree->nodetab[sonj].fathnum == ctree->nodetab[soni].fathnum );
            ctree->nodetab[sonj].fathnum = -2;

            sonsnbr--;
            tmp[j] = tmp[sonsnbr]; j--;
        }
    }
    cnode->sonsnbr = sonsnbr;

    return;
}

/**
 *******************************************************************************
 *
 * @brief Print the elimination tree in a dot file.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 *
 * @param[inout] stream
 *          The file to which write the elimination tree in the dot format.
 *
 *******************************************************************************/
void
candGenDot( const EliminTree *etree,
            const Cand       *candtab,
            FILE             *stream )
{
    pastix_int_t i;

    fprintf(stream,
            "digraph G {\n"
            "\tcolor=white\n"
            "\trankdir=BT;\n");

    for (i=0;  i < etree->nodenbr; i++)
    {
        if ((etree->nodetab[i]).fathnum == -2)
            continue;

        if ( candtab == NULL ) {
            fprintf(stream, "\t\"%ld\" [label=\"#%ld\\nSubtree cost: %e\\nNode cost: %e\\nNode CP: %e\"]\n",
                    (long)i, (long)i,
                    etree->nodetab[i].subtree,
                    etree->nodetab[i].total,
                    etree->nodetab[i].cripath );
        }
        else {
            if ( candtab[i].lcandnum != candtab[i].fcandnum ) {
                fprintf( stream, "\t\"%ld\" [label=\"#%ld\\nCand: %ld - %ld\\nSubtree cost: %e\\nNode cost: %e\\nNode CP: %e\"]\n",
                         (long)i, (long)i,
                         (long)(candtab[i].fcandnum),
                         (long)(candtab[i].lcandnum),
                         etree->nodetab[i].subtree,
                         etree->nodetab[i].total,
                         etree->nodetab[i].cripath );
            }
            else {
                fprintf(stream, "\t\"%ld\" [label=\"#%ld\\nCand: %ld\\nSubtree cost: %e\\nNode cost: %e\\nNode CP: %e\" colorscheme=set312 style=filled fillcolor=%ld]\n",
                        (long)i, (long)i,
                        (long)(candtab[i].fcandnum),
                        etree->nodetab[i].subtree,
                        etree->nodetab[i].total,
                        etree->nodetab[i].cripath,
                        (long)((candtab[i].lcandnum % 12) + 1));
            }
        }
        if ((etree->nodetab[i]).fathnum == -1)
            continue;
        fprintf(stream, "\t\"%ld\"->\"%ld\"\n", (long)i, (long)((etree->nodetab[i]).fathnum));
    }

    fprintf(stream, "}\n");
}

/**
 *******************************************************************************
 *
 * @brief Print the compressed elimination tree in a dot file, where all nodes
 * with the same candidates are merged together.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 *
 * @param[inout] stream
 *          The file to which write the elimination tree in the dot format.
 *
 *******************************************************************************/
void
candGenCompressedDot(const EliminTree *etree, const Cand *candtab, FILE *stream)
{
    EliminTree  *ctree;
    eTreeNode_t *cnode;
    Cand        *ccand;
    pastix_int_t i, cnodesnbr, cnodeidx, *tmp;

    cnodesnbr = compress_getNodeNbr( etree, candtab, eTreeRoot(etree)  );

    fprintf(stderr, "The compressed number of nodes is %ld\n", cnodesnbr);

    /* Let's create a second compressed elimination tree, and the associated candtab */
    MALLOC_INTERN(ctree, 1, EliminTree);
    eTreeInit(ctree);

    MALLOC_INTERN(ccand, cnodesnbr, Cand);
    candInit(ccand, cnodesnbr);

    ctree->nodenbr = cnodesnbr;
    MALLOC_INTERN(ctree->nodetab, cnodesnbr, eTreeNode_t);
    cnode = ctree->nodetab;

    /* Initialize the structure fields */
    for(i=0; i<cnodesnbr; i++, cnode++)
    {
        cnode->total   = 0.;
        cnode->subtree = 0.;
        cnode->cripath = 0.;
        cnode->sonsnbr =  0;
        cnode->fathnum = -1;
        cnode->fsonnum = -1;
    }

    MALLOC_INTERN(tmp, etree->nodenbr, pastix_int_t);
    cnodeidx = -1;
    compress_setSonsNbr( etree, candtab, eTreeRoot( etree ),
                         ctree, ccand, -1, &cnodeidx, tmp );
    memFree_null(tmp);

    /* Write the dot file */
    candGenDot( ctree, ccand, stream );

    /* Free temporary ressources */
    eTreeExit( ctree );
}
