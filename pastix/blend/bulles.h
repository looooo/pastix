#ifndef BULLES_H
#define BULLES_H

/*+ The node structure. +*/
typedef struct BubbleTreeNode_ {
    int                       fathnum;      /*+ index of the father node               +*/
    int                       sonsnbr;
    int                       fsonnum;
    int                       fcandnum;     /*+ first bubble proc                      +*/
    int                       lcandnum;     /*+ last bubble proc                       +*/
    double                    localcost;    /*+ Total cost of the tasks of this node   +*/
    double                    costlevel;    /*+ cost of way from the root to this node +*/
    int                       treelevel;    /*+ cost of way from the root to this node +*/
    int                       nbtasks;      /*+ cost of way from the root to this node +*/
    pastix_int_t              priomin;      /*+ Minimal priority of tasks owned by the bubble +*/
    pastix_int_t              priomax;      /*+ Maximal priority of tasks owned by the bubble +*/
    Queue *                   taskheap;     /*+ Liste de taches de la bulle            +*/
} BubbleTreeNode;


/*+ The bubble tree. +*/
typedef struct BubbleTree_ {
  int                       leavesnbr;    /*+ Number of leaves in tree  +*/
  int                       nodenbr;      /*+ Number of nodes           +*/
  int                       nodemax;      /*+ Number max of nodes       +*/
  int                      *sonstab;      /*+ Array of node             +*/
  BubbleTreeNode           *nodetab;      /*+ Array of node             +*/
} BubbleTree;


#define BFATHER(btree, r)     (btree)->nodetab[r].fathnum
#define BNBSON(btree, r)      (btree)->nodetab[r].sonsnbr
#define BROOT(btree)          (btree)->nodetab[(btree)->leavesnbr]
#define BSON(btree, n, r)     (btree)->sonstab[(btree)->nodetab[n].fsonnum + r ]

void  bubbleInitTree  (BubbleTree *, int);
void  bubbleFree      (BubbleTree *);
int   bubbleAdd       (BubbleTree *, pastix_int_t, pastix_int_t, pastix_int_t, double, double, pastix_int_t);
void  bubbleBuildTree (const BubbleTree *);
void  bubblePrint     (const BubbleTree *, FILE*);

#endif /* BULLES_H */






