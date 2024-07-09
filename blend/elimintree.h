/**
 *
 * @file elimintree.h
 *
 * PaStiX analyse elimin tree and graph header
 *
 * @copyright 1998-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2023-07-21
 *
 * @addtogroup blend_dev_elim
 * @{
 *
 **/
#ifndef _elimintree_h_
#define _elimintree_h_

/**
 * @brief Node of the elimination tree.
 */
typedef struct etree_node_s {
    double       ndecost; /**< Cost of the tree node only (compute + send)   */
    double       ndepath; /**< Critical path of the tree node only           */
    double       subcost; /**< Cost of the subtree (including node)          */
    double       subpath; /**< Critical path of the subtree (including node) */
    int          ndlevel; /**< Node depth in the elimination tree            */
    int          sonsnbr; /**< Number of sons                                */
    pastix_int_t fathnum; /**< index of the father node                      */
    pastix_int_t fsonnum; /**< index of first son                            */
} eTreeNode_t;

/**
 * @brief Elimination tree.
 */
typedef struct etree_s {
    pastix_int_t   baseval; /**< Base value for numberings         */
    pastix_int_t   nodenbr; /**< Number of nodes                   */
    eTreeNode_t  * nodetab; /**< Array of node          [+1,based] */
    pastix_int_t * sonstab; /**< Sons index of nodes               */
} EliminTree;

EliminTree   *eTreeInit      (      pastix_int_t);
void          eTreeExit      (      EliminTree *);
void          eTreeGenDot    (const EliminTree *, FILE *);
void          eTreePrint     (const EliminTree *, FILE *, pastix_int_t );
void          eTreeSetSons   (      EliminTree *);
pastix_int_t  eTreeLeavesNbr (const EliminTree *);
pastix_int_t  eTreeLevel     (const EliminTree *);
pastix_int_t  eTreeNodeLevel (const EliminTree *, pastix_int_t );
EliminTree   *eTreeBuild     (const symbol_matrix_t *);

pastix_int_t eTreeComputeLevels   ( EliminTree *, pastix_int_t, pastix_int_t );
pastix_int_t eTreeGetLevelMinIdx  ( const EliminTree *, pastix_int_t, pastix_int_t, pastix_int_t );

/**
 *******************************************************************************
 *
 * @brief Return the father of a given node.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          Pointer to the elimination tree structure.
 *
 * @param[in] node
 *          The node of interest.
 *
 *******************************************************************************
 *
 * @return The father of the node.
 *
 *******************************************************************************/
static inline eTreeNode_t *
eTreeFather( const EliminTree *etree, pastix_int_t node )
{
    return etree->nodetab + etree->nodetab[node].fathnum;
}

/**
 *******************************************************************************
 *
 * @brief Return the i^{th} son of a given node.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          Pointer to the elimination tree structure.
 *
 * @param[in] node
 *          The node of interest.
 *
 * @param[in] i
 *          The index of the son wanted.
 *
 *******************************************************************************
 *
 * @return The index of the i^{th} son of the node.
 *
 *******************************************************************************/
static inline pastix_int_t
eTreeSonI( const EliminTree *etree, pastix_int_t node, pastix_int_t i )
{
    return etree->sonstab[ etree->nodetab[node].fsonnum + i ];
}

/**
 *******************************************************************************
 *
 * @brief Return the root of the elimination tree.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          Pointer to the elimination tree structure.
 *
 *******************************************************************************
 *
 * @return The index of the root in the elimination tree.
 *
 *******************************************************************************/
static inline pastix_int_t
eTreeRoot( const EliminTree *etree )
{
    (void)etree;
    return -1;
}

#endif /* _elimintree_h_ */

/**
 *@}
 */
