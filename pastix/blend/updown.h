/************************************************************/
/**                                                        **/
/**   NAME       : updown.h                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the UpDown step  .                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     28 oct 1998     **/
/**                                                        **/
/************************************************************/

#ifndef UPDOWN_H
#define UPDOWN_H

/*+ UpDown block structure. +*/

typedef struct UpDownCblk_  {
  PASTIX_INT                       sm2xind;              /*+ Index in the rhs local vector of the unknowns corresponding to the diag blok +*/
  PASTIX_INT *                     browproctab;          /*+ Brow                               +*/
  PASTIX_INT *                     browcblktab;          /*+ Brow                               +*/
  PASTIX_INT                       browprocnbr;          /*+ Brow size                          +*/
  PASTIX_INT                       msgnbr;               /*+ Number of messages                 +*/
  PASTIX_INT volatile              msgcnt;               /*+ Number of messages                 +*/
  PASTIX_INT                       ctrbnbr;              /*+ Number of contributions            +*/
  PASTIX_INT volatile              ctrbcnt;              /*+ Number of contributions            +*/
} UpDownCblk;


/*+ UpDown vector structure. +*/

typedef struct UpDownVector_ {
  UpDownCblk *              cblktab;              /*+ Array of solver column blocks      +*/
  pastix_float_t *                   sm2xtab;              /*+ Unknown vector                     +*/
  PASTIX_INT                       sm2xmax;              /*+ Maximum of coefficients per unknown vector +*/
  PASTIX_INT                       sm2xsze;              /*+ Size of sm2xtab                    +*/
  PASTIX_INT                       sm2xnbr;              /*+ Number of sm2x                     +*/
  PASTIX_INT *                     gcblk2list;           /*+ Global cblknum -> index in listptr +*/
  PASTIX_INT                       gcblk2listnbr;        /*+ Size of gcblk2list                 +*/
  PASTIX_INT *                     listptr;              /*+ Index in list                      +*/
  PASTIX_INT                       listptrnbr;           /*+ Size of listptr                    +*/
  PASTIX_INT *                     listcblk;             /*+ List of cblk in a same row         +*/
  PASTIX_INT *                     listblok;             /*+ List of blok in a same row         +*/
  PASTIX_INT                       listnbr;              /*+ Size of list                       +*/
  PASTIX_INT *                     loc2glob;             /*+ Local cblknum -> global cblknum    +*/
  PASTIX_INT                       loc2globnbr;          /*+ Size of loc2glob                   +*/
  PASTIX_INT *                     lblk2gcblk;           /*+ Local blok -> global facing cblk   +*/
  PASTIX_INT                       gcblknbr;             /*+ total number of cblk               +*/
  PASTIX_INT                       gnodenbr;             /*+ total number of nodes              +*/
  PASTIX_INT                       downmsgnbr;           /*+ Nb messages receive during down    +*/
  PASTIX_INT                       upmsgnbr;             /*+ Nb messages receive during up      +*/
} UpDownVector;

#endif /* UPDOWN_H */
