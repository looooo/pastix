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
  pastix_int_t                       sm2xind;              /*+ Index in the rhs local vector of the unknowns corresponding to the diag blok +*/
  pastix_int_t *                     browproctab;          /*+ Brow                               +*/
  pastix_int_t *                     browcblktab;          /*+ Brow                               +*/
  pastix_int_t                       browprocnbr;          /*+ Brow size                          +*/
  pastix_int_t                       msgnbr;               /*+ Number of messages                 +*/
  pastix_int_t volatile              msgcnt;               /*+ Number of messages                 +*/
  pastix_int_t                       ctrbnbr;              /*+ Number of contributions            +*/
  pastix_int_t volatile              ctrbcnt;              /*+ Number of contributions            +*/
} UpDownCblk;


/*+ UpDown vector structure. +*/

typedef struct UpDownVector_ {
  UpDownCblk *              cblktab;              /*+ Array of solver column blocks      +*/
  pastix_float_t *                   sm2xtab;              /*+ Unknown vector                     +*/
  pastix_int_t                       sm2xmax;              /*+ Maximum of coefficients per unknown vector +*/
  pastix_int_t                       sm2xsze;              /*+ Size of sm2xtab                    +*/
  pastix_int_t                       sm2xnbr;              /*+ Number of sm2x                     +*/
  pastix_int_t *                     gcblk2list;           /*+ Global cblknum -> index in listptr +*/
  pastix_int_t                       gcblk2listnbr;        /*+ Size of gcblk2list                 +*/
  pastix_int_t *                     listptr;              /*+ Index in list                      +*/
  pastix_int_t                       listptrnbr;           /*+ Size of listptr                    +*/
  pastix_int_t *                     listcblk;             /*+ List of cblk in a same row         +*/
  pastix_int_t *                     listblok;             /*+ List of blok in a same row         +*/
  pastix_int_t                       listnbr;              /*+ Size of list                       +*/
  pastix_int_t *                     loc2glob;             /*+ Local cblknum -> global cblknum    +*/
  pastix_int_t                       loc2globnbr;          /*+ Size of loc2glob                   +*/
  pastix_int_t *                     lblk2gcblk;           /*+ Local blok -> global facing cblk   +*/
  pastix_int_t                       gcblknbr;             /*+ total number of cblk               +*/
  pastix_int_t                       gnodenbr;             /*+ total number of nodes              +*/
  pastix_int_t                       downmsgnbr;           /*+ Nb messages receive during down    +*/
  pastix_int_t                       upmsgnbr;             /*+ Nb messages receive during up      +*/
} UpDownVector;

#endif /* UPDOWN_H */
