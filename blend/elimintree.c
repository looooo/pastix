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
 * @copyright 2004-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2024-07-05
 *
 * @addtogroup blend_dev_elim
 * @{
 *
 **/
#include "common.h"
#include "symbol/symbol.h"
#include "blend/elimintree.h"

/**
 *******************************************************************************
 *
 * @brief Initialize the elimination tree structure.
 *
 *******************************************************************************
 *
 * @param[in] nodenbr
 *          TODO
 *
 *******************************************************************************
 *
 * @retval TODO
 *
 *******************************************************************************/
EliminTree *
eTreeInit( pastix_int_t nodenbr )
{
    EliminTree  *etree;
    eTreeNode_t *enode;
    pastix_int_t i;

    MALLOC_INTERN(etree, 1, EliminTree);

    etree->baseval = 0;
    etree->nodenbr = nodenbr;

    MALLOC_INTERN( etree->nodetab, nodenbr + 1, eTreeNode_t);
    MALLOC_INTERN( etree->sonstab, nodenbr,     pastix_int_t);
    memset( etree->sonstab, 0, nodenbr * sizeof(pastix_int_t) );

    /* Initialize the structure fields */
    enode = etree->nodetab;
    for(i=-1; i<nodenbr; i++, enode++)
    {
        enode->ndecost =  0.0;
        enode->ndepath =  0.0;
        enode->subcost =  0.0;
        enode->subpath =  0.0;
        enode->ndlevel = -1;
        enode->sonsnbr =  0;
        enode->fathnum = -1;
        enode->fsonnum = -1;
    }

    /* Shift the nodetab to get the root at -1 position */
    etree->nodetab++;
    return etree;
}

/**
 *******************************************************************************
 *
 * @brief Free the elimination tree structure.
 *
 *******************************************************************************
 *
 * @param[inout] etree
 *          The pointer to the elimination tree to free.
 *
 *******************************************************************************/
void
eTreeExit(EliminTree *etree)
{
    if (etree != NULL) {
        if (etree->nodetab != NULL) {
            etree->nodetab--;
            memFree_null(etree->nodetab);
        }
        memFree_null(etree->sonstab);
        memFree_null(etree);
    }
}

/**
 *******************************************************************************
 *
 * @brief Set the fsonnum fields base on the initialized sonsnbr.
 *
 *******************************************************************************
 *
 * @param[inout] etree
 *          The pointer to the elimination tree for which the fsonnum fields
 *          must be initialized.
 *
 *******************************************************************************/
static inline void
etree_SetSonsIndex(EliminTree *etree)
{
    eTreeNode_t *enode;
    pastix_int_t i;

    /* Set the index of the first sons */
    enode = etree->nodetab - 1;
    enode->fsonnum = 0;
    for(i=0; i<etree->nodenbr; i++, enode++)
    {
        enode[1].fsonnum = enode[0].fsonnum + enode[0].sonsnbr;
    }
    assert((enode[0].fsonnum + enode[0].sonsnbr) == etree->nodenbr);
}

/**
 *******************************************************************************
 *
 * @brief Set the fsonnum fields base on the initialized sonsnbr.
 *
 *******************************************************************************
 *
 * @param[inout] etree
 *          The pointer to the elimination tree for which the fsonnum fields
 *          must be initialized.
 *
 *******************************************************************************/
void
eTreeSetSons(EliminTree *etree)
{
    pastix_int_t i;

    /* Set the index of the first sons */
    etree_SetSonsIndex( etree );

    /* Fill the sonstab */
    for(i=0; i<etree->nodenbr; i++)
    {
        eTreeNode_t *efather = eTreeFather(etree, i);
        pastix_int_t node = efather->fsonnum;
        efather->fsonnum++;
        assert( (node >= 0) && (node < etree->nodenbr) );
        etree->sonstab[ node ] = i;
    }

    /* Restore fsonnum fields */
    etree_SetSonsIndex( etree );
}

/**
 *******************************************************************************
 *
 * @brief Compute the number of leaves.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 *******************************************************************************
 *
 * @return The number of leaves in the elimination tree.
 *
 *******************************************************************************/
pastix_int_t
eTreeLeavesNbr(const EliminTree *etree)
{
    pastix_int_t i;
    pastix_int_t leavenbr;
    leavenbr = 0;
    for(i=0;i<etree->nodenbr;i++) {
        if(etree->nodetab[i].sonsnbr == 0) {
            leavenbr++;
        }
    }

    return leavenbr;
}

/**
 *******************************************************************************
 *
 * @brief Compute the height of the elimination tree.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 *******************************************************************************
 *
 * @return The height of the elimination tree.
 *
 *******************************************************************************/
pastix_int_t
eTreeLevel(const EliminTree *etree)
{
    pastix_int_t maxlevel;
    pastix_int_t nodelevel;
    pastix_int_t i;
    maxlevel = 0;
    for(i=0;i<etree->nodenbr;i++)
    {
        nodelevel = eTreeNodeLevel(etree, i);
        if(nodelevel > maxlevel) {
            maxlevel = nodelevel;
        }
    }

    return maxlevel;
}

/**
 *******************************************************************************
 *
 * @brief Compute the number of level existing below a given node.
 *
 *******************************************************************************
 *
 * @param[inout] etree
 *          The pointer to the elimination tree. On exit, all traversed nodes
 *          will be updated with the depth related to the given root.
 *
 *
 * @param[in] rootnum
 *          The root of the tree to study.
 *
 * @param[in] rootlvl
 *          The level of the root node.
 *
 *******************************************************************************
 *
 * @return The maximal depth of the tree.
 *
 *******************************************************************************/
pastix_int_t
eTreeComputeLevels( EliminTree *etree, pastix_int_t rootnum, pastix_int_t rootlvl )
{
    pastix_int_t sonsnbr, i, son, max;
    pastix_int_t maxlevel;

    etree->nodetab[ rootnum ].ndlevel = rootlvl;
    sonsnbr = etree->nodetab[ rootnum ].sonsnbr;

    maxlevel = rootlvl;
    rootlvl++;
    for(i=0;i<sonsnbr;i++)
    {
        son = eTreeSonI( etree, rootnum, i );
        max = eTreeComputeLevels( etree, son, rootlvl );
        maxlevel = pastix_imax( max, maxlevel );
    }
    return maxlevel;
}

/**
 *******************************************************************************
 *
 * @brief Compute the number of level existing below a given node.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 *
 * @param[in] nodenum
 *          The index of the node to study.
 *
 *******************************************************************************
 *
 * @return The number of level below the node including it.
 *
 *******************************************************************************/
pastix_int_t
eTreeNodeLevel(const EliminTree *etree, pastix_int_t nodenum )
{
    pastix_int_t level;

    level = 0;
    /* If fake root no levels */
    if(nodenum == eTreeRoot(etree)) {
        return 0;
    }
    level++;
    while(etree->nodetab[nodenum].fathnum != eTreeRoot(etree))
    {
        level++;
        nodenum = etree->nodetab[nodenum].fathnum;
    }
    return level;
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
eTreeGenDot( const EliminTree *etree,
             FILE             *stream )
{
    pastix_int_t i;

    fprintf(stream,
            "digraph G {\n"
            "\tcolor=white\n"
            "rankdir=BT;\n");

    for (i=0;  i < etree->nodenbr; i++)
    {
        fprintf(stream, "\t\"%ld\" [label=\"#%ld\\nSubtree cost: %e\\nNode cost: %e\\nNode CP: %e\"]\n",
                (long)i, (long)i,
                etree->nodetab[i].subcost,
                etree->nodetab[i].ndepath,
                etree->nodetab[i].subpath );

        if ((etree->nodetab[i]).fathnum == -1) {
            continue;
        }
        fprintf( stream, "\t\"%ld\"->\"%ld\"\n",
                 (long)i,
                 (long)((etree->nodetab[i]).fathnum) );
    }

    fprintf(stream, "}\n");
}

/**
 *******************************************************************************
 *
 * @brief Print the elimination tree in a human readable format.
 *
 * Each node is writen as:
 *  Rootnum idx number_of_sons:
 *          (son_1)
 *          (son_2)
 *          ...
 *          (son_n)
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
 * @param[in] rootnum
 *          The root of the subtree to write into the file.
 *
 *******************************************************************************/
void
eTreePrint(const EliminTree *etree, FILE *stream, pastix_int_t rootnum )
{
    int i, sonsnbr;
    pastix_int_t son;

    sonsnbr = etree->nodetab[ rootnum ].sonsnbr;

    fprintf(stream, "Rootnum %ld %d\n", (long)rootnum, sonsnbr);
    for(i=0;i<sonsnbr;i++) {
        fprintf(stream,"       (%4ld)\n",  (long)eTreeSonI(etree, rootnum, i));
    }

    for(i=0;i<sonsnbr;i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        if (etree->nodetab[son].sonsnbr) {
            eTreePrint(etree, stream, son);
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Return the smallest index of the nodes belonging to the level lvl in
 * the elimination tree.
 *
 * Look recursively in the elimination tree to return the smallest index of all
 * the nodes belonging to a same level defined by lvl.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 * @param[in] root
 *          The root of the subtree to study.
 *
 * @param[in] lvl
 *          The number of levels below root to look for.
 *
 * @param[in] idxmin
 *          The minimum index already discovered.
 *
 *******************************************************************************
 *
 * @return The minimum index found in the descendants of the root node, lvl
 * levels below.
 *
 *******************************************************************************/
pastix_int_t
eTreeGetLevelMinIdx( const EliminTree *etree,
                     pastix_int_t      root,
                     pastix_int_t      lvl,
                     pastix_int_t      idxmin )
{
    int i, sonsnbr;
    pastix_int_t idx, son;

    idx = (root == -1) ? idxmin : pastix_imin( root, idxmin );

    if ( lvl == 0 ) {
        return idx;
    }

    sonsnbr = etree->nodetab[ root ].sonsnbr;
    if ( sonsnbr == 0 ) {
        return idx;
    }

    lvl--;
    for(i=0;i<sonsnbr;i++)
    {
        son = eTreeSonI(etree, root, i);
        idx = eTreeGetLevelMinIdx( etree, son, lvl, idx );
    }
    return idx;
}

/**
 *******************************************************************************
 *
 * @brief Build the elimination tree.
 *
 * The elimination tree is computed based on a given symbolic structure, and
 * not from the tree given by the ordering library. Each father of a node is
 * defined as the facing column block of the first off diagonal block.
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          The pointer to the symbol matrix from which the elimination tree is
 *          computed.
 *
 *******************************************************************************
 *
 * @return The elimination tree linked to the symbol matrix.
 *
 *******************************************************************************/
EliminTree *
eTreeBuild(const symbol_matrix_t *symbmtx)
{
    eTreeNode_t *enode, *eroot;
    EliminTree *etree = NULL;
    pastix_int_t i;
    pastix_int_t totalsonsnbr;

    etree = eTreeInit( symbmtx->cblknbr );
    eroot = &(etree->nodetab[eTreeRoot(etree)]);

    /* Compute the fathers and the number of sons */
    totalsonsnbr = 0;
    enode = etree->nodetab;
    for(i=0; i<symbmtx->cblknbr; i++, enode++)
    {
        /* If the cblk has at least one extra diagonal block,          */
        /* the father of the node is the facing block of the first odb */
        if( (symbmtx->cblktab[i+1].bloknum - symbmtx->cblktab[i].bloknum) > 1 )
        {
            pastix_int_t fathnum = symbmtx->bloktab[ symbmtx->cblktab[i].bloknum+1 ].fcblknm;
            enode->fathnum = fathnum;
            etree->nodetab[fathnum].sonsnbr++;
            totalsonsnbr++;
        }
        /* Otherwise this is a root, and we attach it to -1 */
        else {
            etree->nodetab[i].fathnum = -1;
            eroot->sonsnbr++;
            totalsonsnbr++;
        }
    }

    /* Check that we have only one root */
    assert(totalsonsnbr == symbmtx->cblknbr);

    /* Set the sons */
    eTreeSetSons( etree );

    return etree;
}

/**
 *@}
 */
