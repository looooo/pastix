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

/**
 *******************************************************************************
 *
 * @brief Initialize the elimination tree structure.
 *
 *******************************************************************************
 *
 * @param[inout] etree
 *          The pointer to the allocated structure to initialize.
 *
 *******************************************************************************/
void
eTreeInit(EliminTree *etree)
{
    etree->baseval = 0;
    etree->nodenbr = 0;
    etree->nodetab = NULL;
    etree->sonstab = NULL;
    return;
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
    memFree_null(etree->nodetab);
    memFree_null(etree->sonstab);
    memFree_null(etree);
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
    for(i=0;i<etree->nodenbr;i++)
        if(etree->nodetab[i].sonsnbr == 0)
            leavenbr++;

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
        if(nodelevel>maxlevel)
            maxlevel = nodelevel;
    }

    return maxlevel;
}

/**
 *******************************************************************************
 *
 * @brief Compute the nuber of level existing below a given node.
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

    level = 1;
    if(nodenum == eTreeRoot(etree))
        return level;
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
 * @param[in] out
 *          The file to which write the elimination tree in the dot format.
 *
 *******************************************************************************/
void
eTreeGenDot(const EliminTree *etree, FILE *out)
{
    pastix_int_t i;

    fprintf(out,
            "digraph G {\n"
            "\tcolor=white\n"
            "rankdir=BT;\n");

    for (i=0;  i < etree->nodenbr; i++)
    {
        if ((etree->nodetab[i]).fathnum == -1)
            continue;
        fprintf(out, "\t\"%ld\"->\"%ld\"\n", (long)i, (long)((etree->nodetab[i]).fathnum));
    }

    fprintf(out, "}\n");
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
 * @param[in] out
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
    for(i=0;i<sonsnbr;i++)
        fprintf(stream,"       (%4ld)\n",  (long)eTreeSonI(etree, rootnum, i));

    for(i=0;i<sonsnbr;i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        if (etree->nodetab[son].sonsnbr)
            eTreePrint(etree, stream, son);
    }
}

/**
 *******************************************************************************
 *
 * @brief Build the elimination tree.
 *
 * The elimination tree is computed based on a given symbvolic structure,a nd
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
eTreeBuild(const SymbolMatrix *symbmtx)
{
    eTreeNode_t *enode;
    EliminTree *etree = NULL;
    pastix_int_t i;
    pastix_int_t totalsonsnbr;
    pastix_int_t sonstabcur;

    MALLOC_INTERN(etree, 1, EliminTree);
    eTreeInit(etree);

    etree->nodenbr = symbmtx->cblknbr;
    MALLOC_INTERN(etree->nodetab, etree->nodenbr, eTreeNode_t);
    enode = etree->nodetab;

    /* Initialize the structure fields */
    for(i=0; i<symbmtx->cblknbr; i++, enode++)
    {
        enode->total   = 0.;
        enode->subtree = 0.;
        enode->sonsnbr =  0;
        enode->fathnum = -1;
        enode->fsonnum = -1;
    }

    totalsonsnbr = 0;
    for(i=0; i<symbmtx->cblknbr; i++)
    {
        /* If the cblk has at least one extra diagonal block,          */
        /* the father of the node is the facing block of the first odb */
        if( (symbmtx->cblktab[i+1].bloknum - symbmtx->cblktab[i].bloknum) > 1 )
        {
            etree->nodetab[i].fathnum = symbmtx->bloktab[ symbmtx->cblktab[i].bloknum+1 ].fcblknm;
            (eTreeFather( etree, i )->sonsnbr)++;
            totalsonsnbr++;
        }
#if defined(PASTIX_DEBUG_BLEND)
        else
        {
            if(i != (symbmtx->cblknbr-1)) {
                fprintf(stderr, "Cblk %ld has no extradiagonal %ld %ld !! \n", (long)i,
                        (long)symbmtx->cblktab[i].bloknum, (long)symbmtx->cblktab[i+1].bloknum);
                assert( 0 );
            }
        }
#endif
    }

    /* Check that we have only one root */
    assert(totalsonsnbr == (symbmtx->cblknbr-1));

    if( totalsonsnbr > 0 ) {
        MALLOC_INTERN(etree->sonstab, totalsonsnbr, pastix_int_t);
    }

    /* Set the index of the first sons */
    sonstabcur = 0;
    for(i=0; i<symbmtx->cblknbr; i++)
    {
        etree->nodetab[i].fsonnum = sonstabcur;
        sonstabcur += etree->nodetab[i].sonsnbr;
    }
    assert(sonstabcur == totalsonsnbr);

    /* Fill the sonstab */
    /* No need to go to the root */
    for(i=0; i<symbmtx->cblknbr-1; i++)
    {
        etree->sonstab[ (eTreeFather(etree, i)->fsonnum)++] = i;
    }

    /* Restore fsonnum fields */
    sonstabcur = 0;
    for(i=0; i<symbmtx->cblknbr; i++)
    {
        etree->nodetab[i].fsonnum = sonstabcur;
        sonstabcur += etree->nodetab[i].sonsnbr;
    }
    assert(sonstabcur == totalsonsnbr);

    return etree;
}

/**
 *@}
 */
