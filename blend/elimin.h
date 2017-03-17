/**
 *
 * @file elimin.h
 *
 * PaStiX analyse elimin tree and graph header
 *
 * @copyright 1998-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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

/**
 * @brief Node of the elimination tree.
 */
typedef struct eTreeNode_s {
    double         total;   /**< Cost of the treenode only (compute + send) */
    double         subtree; /**< Cost of the subtree (includes total)       */
    pastix_int_t   sonsnbr; /**< Number of sons                             */
    pastix_int_t   fathnum; /**< index of the father node                   */
    pastix_int_t   fsonnum; /**< index of first son                         */
} eTreeNode_t;

/**
 * @brief Elimination tree.
 */
typedef struct EliminTree_ {
    pastix_int_t   baseval; /**< Base value for numberings         */
    pastix_int_t   nodenbr; /**< Number of nodes                   */
    eTreeNode_t  * nodetab; /**< Array of node          [+1,based] */
    pastix_int_t * sonstab; /**< Sons index of nodes               */
} EliminTree;

pastix_int_t  eTreeInit      (      EliminTree *);
void          eTreeExit      (      EliminTree *);
void          eTreeGenDot    (const EliminTree *, FILE *);
void          eTreePrint     (const EliminTree *, FILE *, pastix_int_t );
pastix_int_t  eTreeLeavesNbr (const EliminTree *);
pastix_int_t  eTreeLevel     (const EliminTree *);
pastix_int_t  eTreeNodeLevel (const EliminTree *, pastix_int_t );
EliminTree   *eTreeBuild     (const SymbolMatrix *);

#define eTreeFather( __etree__, __node__ )       ((__etree__)->nodetab[(__etree__)->nodetab[(__node__)].fathnum])
#define eTreeSonI( __etree__, __node__, __i__ )  ((__etree__)->sonstab[(__etree__)->nodetab[(__node__)].fsonnum + (__i__)])
#define eTreeRoot( __etree__ )                   ((__etree__)->nodenbr - 1)
