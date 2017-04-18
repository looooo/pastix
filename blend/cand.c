/**
 *
 * @file cand.c
 *
 * PaStiX analyse functions to manipulate candidates on the elimination tree
 * structure.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @addtogroup blend_dev_cand
 * @{
 *    This module contains all subroutines to initialize the candidates array
 *    for each supernode, as well as supernode properties that are defined by
 *    level such as 2D layouts and 2D tasks.
 *
 **/
#include "common.h"
#include "symbol.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "solver.h"

/**
 *******************************************************************************
 *
 * @brief Initialize the candtab array with default values.
 *
 *******************************************************************************
 *
 * @param[inout] candtab
 *          The array of size cblknbr of Cand structure to initialize.
 *
 * @param[in] cblknbr
 *          The size of the candtab array.
 *
 *******************************************************************************/
void
candInit( Cand *candtab,
          pastix_int_t cblknbr )
{
    pastix_int_t i;
    for(i=0;i<cblknbr;i++)
    {
        candtab[i].costlevel = 0.0;
        candtab[i].treelevel = 0;
        candtab[i].fcandnum  = -1;
        candtab[i].lcandnum  = -1;
        candtab[i].fccandnum = -1;
        candtab[i].lccandnum = -1;
        candtab[i].cluster   = -1;
        candtab[i].cblktype  = CBLK_LAYOUT_2D | CBLK_TASKS_2D;
    }
}

/**
 *******************************************************************************
 *
 * @brief Print the candidates array into the candtab.txt file
 *
 *******************************************************************************
 *
 * @param[in] candtab
 *          The array of size cblknbr to print in the file.
 *
 * @param[in] cblknbr
 *          The size of the candtab array.
 *
 *******************************************************************************/
void
candSave( const Cand *candtab,
          pastix_int_t cblknbr )
{
    FILE *f = fopen("candtab.txt", "w");
    pastix_int_t i;
    fprintf(f, "%ld\n", cblknbr );
    for(i=0;i<cblknbr;i++)
    {
        fprintf(f, "%lf %ld %ld %ld %ld %ld %ld %d\n",
                candtab[i].costlevel,
                candtab[i].treelevel,
                candtab[i].fcandnum,
                candtab[i].lcandnum,
                candtab[i].fccandnum,
                candtab[i].lccandnum,
                candtab[i].cluster,
                candtab[i].cblktype );
    }
    fclose(f);
}

/**
 *******************************************************************************
 *
 * @brief Set a single candidate recursively to a subtree
 *
 *******************************************************************************
 *
 * @param[inout] candtab
 *          On entry, the array of candidates initialized with candInit().
 *          On exit, each node belonging to the subtree has procnum as a single
 *          candidate.
 *
 * @param[in] etree
 *          The full elimination tree structure of the problem
 *
 * @param[in] rootnum
 *          The root index of the subtree to initialize
 *
 * @param[in] procnum
 *          The processor index to affect to all nodes in the subtree
 *
 *******************************************************************************/
void
candSetSubCandidate( Cand *candtab,
                     const EliminTree *etree,
                     pastix_int_t rootnum,
                     pastix_int_t procnum )
{
    pastix_int_t i;

    candtab[rootnum].fcandnum = procnum;
    candtab[rootnum].lcandnum = procnum;

    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
        candSetSubCandidate( candtab, etree, eTreeSonI(etree, rootnum, i), procnum );
}

/**
 *******************************************************************************
 *
 * @brief Set the clusters candidates from the cores canditates
 *
 *******************************************************************************
 *
 * @param[inout] candtab
 *          On entry, the array of candidates with the first and last core
 *          candidates initialized.
 *          On exit, the array of candidate with te first and last cluster
 *          candidate information updated.
 *
 * @param[in] cblknbr
 *          The size of the candtab array.
 *
 * @param[in] core2clust
 *          An array that defines the cluster (MPI process) that owns each core
 *          candidate.
 *
 * @param[in] coresnbr
 *          The size of the core2clust array.
 *
 *******************************************************************************/
void
candSetClusterCand(       Cand         *candtab,
                          pastix_int_t  cblknbr,
                    const pastix_int_t *core2clust,
                          pastix_int_t  coresnbr )
{
    pastix_int_t i;
    (void)coresnbr;

    for(i=0; i<cblknbr; i++) {
        assert( candtab[i].fcandnum >= 0 );
        assert( candtab[i].lcandnum >= 0 );
        assert( candtab[i].fcandnum < coresnbr );
        assert( candtab[i].lcandnum < coresnbr );
        candtab[i].fccandnum = core2clust[ candtab[i].fcandnum ];
        candtab[i].lccandnum = core2clust[ candtab[i].lcandnum ];
    }
}

/**
 *******************************************************************************
 *
 * @brief Check the correctness of the computed candidates
 *
 * Each node of the elimination tree must have a set of candidates included in
 * its father's set.
 *
 *******************************************************************************
 *
 * @param[in] candtab
 *          On entry, the array of candidates to check.
 *
 * @param[in] symbmtx
 *          The symbol matrix structure associated to the candidate array.
 *
 *******************************************************************************
 *
 * @retval 0 if bad candidat set appear.
 * @retval 1 if success.
 *
 *******************************************************************************/
int
candCheck( const Cand         *candtab,
           const SymbolMatrix *symbmtx )
{
    pastix_int_t i, j;
    pastix_int_t facecblknum;

    for(i=0; i<symbmtx->cblknbr; i++)
    {
        for(j = symbmtx->cblktab[i].bloknum;
            j < symbmtx->cblktab[i+1].bloknum; j++)
        {
            facecblknum = symbmtx->bloktab[j].fcblknm;

            if( (candtab[i].fcandnum < candtab[facecblknum].fcandnum) ||
                (candtab[i].lcandnum > candtab[facecblknum].lcandnum) )
            {
                errorPrint("bad processor candidat sets : cblk %ld candidat =[%ld %ld] father %ld candidat = [%ld %ld].",
                           (long)i, (long)candtab[i].fcandnum, (long)candtab[i].lcandnum,
                           (long)facecblknum, (long)candtab[facecblknum].fcandnum,
                           (long)candtab[facecblknum].lcandnum);
                return 0;
            }
        }
    }
    return 1;
}

/**
 *******************************************************************************
 *
 * @brief Recursive function to update the cost fields of the both the candtab
 *        array, and the elimination tree structure.
 *
 *******************************************************************************
 *
 * @param[in]    rootnum
 *               Root of the subtree.
 *
 * @param[inout] candtab
 *               Pointer to the global candtab array where fields treelevel and
 *               costlevel are updated. Treelevel represent the depth of the
 *               node in the elimination tree, and costlevel the cost of the
 *               path from rootnum to the root of the tree.
 *
 * @param[inout] etree
 *               Pointer to the global elimination tree structure. The node
 *               fields total and subtree are updated. Total represents the
 *               total cost of of the current node only, and subtree the cost of
 *               the current node and all its descendents.
 *
 * @param[in]    symbmtx
 *               Pointer to the symbol matrix we are working with.
 *
 * @param[in]    costmtx
 *               Pointer to the cost matrix associated to the symbol matrix and
 *               that holds the cost of each cblk and blok.
 *
 *******************************************************************************
 *
 * @return The cost of the subtree.
 *
 *******************************************************************************/
static inline double
candSubTreeBuild( pastix_int_t        rootnum,
                  Cand               *candtab,
                  EliminTree         *etree,
                  const SymbolMatrix *symbmtx,
                  const CostMatrix   *costmtx )
{
    double cost;
    pastix_int_t i, son, bloknum;

    /* Compute cost of current node */
#if defined(BLEND_COST_LL)
    cost = costmtx->cblkcost[rootnum];
#else
    cost = 0.;
    for( bloknum = symbmtx->cblktab[ rootnum   ].bloknum;
         bloknum < symbmtx->cblktab[ rootnum+1 ].bloknum; bloknum++)
    {
        cost += costmtx->blokcost[ bloknum ];
    }
#endif

    etree->nodetab[ rootnum ].total   = cost;
    etree->nodetab[ rootnum ].subtree = cost;

    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        candtab[ son ].treelevel = candtab[ rootnum ].treelevel - 1;
        candtab[ son ].costlevel = candtab[ rootnum ].costlevel - cost;

        etree->nodetab[ rootnum ].subtree +=
            candSubTreeBuild( son, candtab, etree, symbmtx, costmtx );
    }

    return etree->nodetab[ rootnum ].subtree;
}

/**
 *******************************************************************************
 *
 * @brief Recursive function to compute the distribution of the nodes among the
 * different levels.
 *
 * This function defines which cblk are candidates to be stored with a 2D layout,
 * computed as 2D tasks, and/or be part of the Schur complement. The criteria to
 * remove the 2D task flags is the width of the cblk.
 *
 *******************************************************************************
 *
 * @param[in]    rootnum
 *               Root of the subtree.
 *
 * @param[in]    cblktype
 *               List of flags that can be forwarded to rootnum and its
 *               descendents.
 *
 * @param[in]    ratiolimit
 *               Ratio that defines the minimal size to allow the flag 2D tasks
 *               to be forwarded to the sons.
 *
 * @param[inout] candtab
 *               Pointer to the global candtab array where field cblktype is
 *               updated. cblktype defines the optimization/properties that are
 *               defined on each cblk and which are defined by level in the
 *               tree.
 *
 * @param[in]    etree
 *               Pointer to the global elimination tree structure that is used
 *               to travel through the cblk, and affect the properies with the
 *               correct filiation property.
 *
 * @param[in]    symbmtx
 *               Pointer to the symbol matrix we are working with.
 *
 * @param[in]    costmtx
 *               Pointer to the cost matrix associated to the symbol matrix and
 *               that holds the cost of each cblk and blok.
 *
 *******************************************************************************/
static inline void
candSubTreeDistribWithSize( pastix_int_t        rootnum,
                            pastix_int_t        cblktype,
                            pastix_int_t        ratiolimit,
                            Cand               *candtab,
                            const EliminTree   *etree,
                            const SymbolMatrix *symbmtx )
{
    pastix_int_t i, son;

    if ( (cblktype & CBLK_IN_SCHUR) &&
         (symbmtx->cblktab[ rootnum ].lcolnum < symbmtx->schurfcol) )
    {
        cblktype &= ~(CBLK_IN_SCHUR);
    }

    if(cblktype & CBLK_TASKS_2D) {
        pastix_int_t width = symbmtx->cblktab[ rootnum ].lcolnum - symbmtx->cblktab[ rootnum ].fcolnum + 1;

        if((ratiolimit >= 0) && (width >= ratiolimit))
            candtab[ rootnum ].cblktype = cblktype;
        else
            candtab[ rootnum ].cblktype = cblktype & (~CBLK_TASKS_2D);
    }
    else {
        candtab[ rootnum ].cblktype = cblktype;
    }

    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        candSubTreeDistribWithSize( son, candtab[ rootnum ].cblktype, ratiolimit, candtab, etree, symbmtx);
    }
}

/**
 *******************************************************************************
 *
 * @brief Recursive function to compute the distribution of the nodes among the
 * different levels based on depth.
 *
 * This function defines which cblk are candidates to be stored with a 2D layout,
 * computed as 2D tasks, and/or be part of the Schur complement. The criteria to
 * remove the 2D task flags is the depth on the elimination tree.
 *
 *******************************************************************************
 *
 * @param[in]    depth
 *               Depth until which 2D tasks is enabled.
 *
 * @param[in]    symbmtx
 *               Pointer to the symbol matrix we are working with.
 *
 * @param[inout] candtab
 *               Pointer to the global candtab array where the field cblktype is
 *               updated. cblktype defines the optimization/properties that are
 *               defined on each cblk and which are defined by level in the
 *               tree.
 *
 *******************************************************************************/
static inline void
candDistribWithDepth( pastix_int_t        depth,
                      const SymbolMatrix *symbmtx,
                      Cand               *candtab )
{
    pastix_int_t i, cblknbr, schurfcol;

    cblknbr   = symbmtx->cblknbr;
    schurfcol = symbmtx->schurfcol;

    for(i=0;i<cblknbr;i++)
    {
        if( candtab[i].treelevel > depth )
            candtab[i].cblktype = CBLK_LAYOUT_2D | CBLK_TASKS_2D;
        else
            candtab[i].cblktype = CBLK_LAYOUT_2D;

        if ( symbmtx->cblktab[ i ].fcolnum >= schurfcol ) {
            candtab[i].cblktype &= CBLK_IN_SCHUR;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Finish to build the candtab array for the proportionnal mapping.
 *
 * This function defines which cblk are candidates to be stored with a 2D
 * layout, computed as 2D tasks, and/or be part of the Schur complement. It also
 * set the cost of each cblk and the total cost of each subtree.
 *
 *******************************************************************************
 *
 * @param[in]    autolevel
 *               If autolevel is defined, the width criteria is used to defined
 *               the 2D tasks flag on the cblk.
 *
 * @param[in]    level2D
 *               If autolevel is false, then level2d defines te number of levels
 *               handles as 2D tasks.
 *
 * @param[in]    ratiolimit
 *               if autolevel is true, ratiolimit defines the minimal width of
 *               the node that will b e handled as 2D tasks.
 *
 * @param[inout] candtab
 *               Pointer to the global candtab array that needs to be initialized.
 *
 * @param[in]    symbmtx
 *               Pointer to the symbol matrix we are working with.
 *
 * @param[in]    costmtx
 *               Pointer to the cost matrix associated to the symbol matrix and
 *               that holds the cost of each cblk and blok.
 *
 *******************************************************************************/
void
candBuild( pastix_int_t autolevel, pastix_int_t level2D, pastix_int_t ratiolimit,
           Cand               *candtab,
           EliminTree         *etree,
           const SymbolMatrix *symbmtx,
           const CostMatrix   *costmtx )
{
    pastix_int_t root = eTreeRoot(etree);

    /* Let's start with the root */
    candtab[ root ].costlevel = -1.0;
    candtab[ root ].treelevel = -1;

    candSubTreeBuild( root, candtab, etree, symbmtx, costmtx );

    /* Let's set the cblk type of each node */
    if(autolevel)
    {
        candSubTreeDistribWithSize( eTreeRoot(etree),
                                    CBLK_LAYOUT_2D | CBLK_TASKS_2D | CBLK_IN_SCHUR,
                                    ratiolimit,
                                    candtab, etree, symbmtx );
    }
    else
    {
        candDistribWithDepth( -level2D, symbmtx, candtab );
    }
}

/**
 *@}
 */

