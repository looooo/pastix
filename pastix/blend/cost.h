/*+********************************************************+*/
/*+                                                        +*/
/*+   NAME       : cost.h                                  +*/
/*+                                                        +*/
/*+   AUTHORS    : Pascal HENON                            +*/
/*+                                                        +*/
/*+   FUNCTION   : Part of a parallel direct block solver. +*/
/*+                Definitions of the structures that      +*/
/*+                contains cost informations              +*/
/*+                                                        +*/
/*+   DATES      : # Version 0.0  : from : 27 sep 1998     +*/
/*+                                 to     03 sep 1998     +*/
/*+                                                        +*/
/*+********************************************************+*/

/*
**  The type and structure definitions.
*/

typedef struct CostCblk_ {
  double                     send;    /* Communication cost                         */
  double                     compute; /* Compute cost                               */
  double                     total;   /* Cost of the treenode only (compute + send) */
  double                     subtree; /* Cost of the subtree (included total)       */

} CostCblk;

typedef struct CostBlok_ {
  double                    contrib; /*+ Cost of contrib bounded to this blok                  +*/
  PASTIX_INT                       linenbr; /*+ Number of no empty line above the blok (blok include) +*/
} CostBlok;

typedef struct CostMatrix_ {
    CostCblk              *     cblktab;
    CostBlok              *     bloktab;
} CostMatrix;

/*
**  The function prototypes.
*/

#ifndef COST
#define static
#endif
PASTIX_INT                         costInit            (CostMatrix *);
void                        costExit            (CostMatrix *);
PASTIX_INT                         costLoad            (CostMatrix *, FILE *);
PASTIX_INT                         costSave            (CostMatrix *, FILE *);

#undef static
