/************************************************************/
/**                                                        **/
/**   NAME       : elimin.h                                **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the elimination tree.               **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     27 jul 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ The node structure. +*/
typedef struct TreeNode_ {
  pastix_int_t                       sonsnbr;              /*+ Number of sons                          +*/
  pastix_int_t                       fathnum;              /*+ index of the father node                +*/
  pastix_int_t                       fsonnum;              /*+ index of first son                      +*/
} TreeNode;

/*+ The elimination tree. +*/

typedef struct EliminTree_ {
  pastix_int_t                       baseval;              /*+ Base value for numberings         +*/
  pastix_int_t                       nodenbr;              /*+ Number of nodes                   +*/
  TreeNode   *              nodetab;              /*+ Array of node          [+1,based] +*/
  pastix_int_t        *              sonstab;              /*+ Sons index of nodes               +*/  
} EliminTree;


/*+ The elimination graph. +*/
/*+ we only need the in-edges between graph vertex
    Out-edges can be found with the symbol matrix data +*/
/* OIMBE innbr ne sert pas necessairement !*/
typedef struct EliminVertex_ {
  pastix_int_t                       innum;                /*+ index of first in-bloc            +*/
  pastix_int_t                       innbr;                /*+ number of in-blocs                +*/   
} EliminVertex;

typedef struct EliminGraph_ {
  pastix_int_t                       baseval;              /*+ Base value for numberings         +*/
  pastix_int_t                       vertnbr;              /*+ number of vertex in the graph     +*/
  EliminVertex   *          verttab;              /*+ Array of vertex                   +*/           
  pastix_int_t            *          inbltab;              /*+ Array of in-blocs index           +*/
  pastix_int_t            *          ownetab;              /*+ Array of cbl owner bloc           +*/           
} EliminGraph;

/*
**  The function prototypes.
*/

#ifndef STRUCT_ELIMINTREE
#define static
#endif
pastix_int_t                         egraphInit        (EliminGraph *);
void                        egraphExit        (EliminGraph *);
pastix_int_t                         egraphLoad        (EliminGraph *, FILE *);
pastix_int_t                         egraphSave        (EliminGraph *, FILE *);
pastix_int_t                         treeInit          (EliminTree *);
void                        treeExit          (EliminTree *);
pastix_int_t                         treeLoad          (EliminTree *, FILE *);
pastix_int_t                         treeSave          (EliminTree *, FILE *);
void                        treePlot          (EliminTree *, FILE *);

#undef static

#define TFATHER(treenode, r)     treenode->nodetab[treenode->nodetab[r].fathnum]
/* return the i_th sons num of the node n */ 
#define TSON(treenode, n, r)     treenode->sonstab[ treenode->nodetab[n].fsonnum + r ]
#define ROOT(treenode)           treenode->nodenbr-1
