/************************************************************/
/**                                                        **/
/**   NAME       : extrastruct.h                           **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Structures used to add extra cblk       **/
/**                when a cblk is splitted.                **/
/**                The goal is to avoid a big amount of    **/
/**                dynamic memory allocations and          **/
/**                memory copies                           **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     08 sep 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/
typedef struct ExtraCostMatrix_ {
  CostCblk   *              cblktab;
  CostBlok   *              bloktab;
} ExtraCostMatrix;



typedef struct ExtraSymbolMatrix_ {
  PASTIX_INT                       baseval;              /*+ Base value for numberings                    +*/
  PASTIX_INT                       addcblk;              /*+ Number of cblk created                       +*/
  PASTIX_INT                       addblok;              /*+ Number of blok created                       +*/
  PASTIX_INT        *              sptcblk;              /*+ Index for splitted cblk in the cblktab       +*/
  PASTIX_INT        *              sptcbnb;              /*+ Number of splitted cblk for a cblk           +*/
  PASTIX_INT        *              sptblok;              /*+ Index for splitted blok in the bloktab       +*/
  PASTIX_INT        *              sptblnb;              /*+ Number of splitted blok for a blok           +*/
  PASTIX_INT        *              subtreeblnbr;         /*+ Number of blok in the subtree                +*/
  PASTIX_INT                       curcblk;              /*+ Cursor for cblktab                           +*/
  PASTIX_INT                       sizcblk;              /*+ Size of allocated cblktab                    +*/
  PASTIX_INT                       curblok;              /*+ Cursor for bloktab                           +*/
  PASTIX_INT                       sizblok;              /*+ Size of allocated bloktab                    +*/
  SymbolCblk *              cblktab;              /*+ Array of column blocks [+1,based]            +*/
  SymbolBlok *              bloktab;              /*+ Array of blocks [based]                      +*/
} ExtraSymbolMatrix;

/*
**  The function prototypes.
*/

#ifndef EXTRASYMBOL
#define static
#endif
void                        extra_inc_cblk           (ExtraSymbolMatrix *, ExtraCostMatrix *);
void                        extra_inc_blok           (ExtraSymbolMatrix *, ExtraCostMatrix *);
PASTIX_INT                         extrasymbolInit          (ExtraSymbolMatrix *);
void                        extrasymbolExit          (ExtraSymbolMatrix *);
PASTIX_INT                         extrasymbolLoad          (ExtraSymbolMatrix *, FILE *);
PASTIX_INT                         extrasymbolSave          (ExtraSymbolMatrix *, FILE *);
PASTIX_INT                         extracostInit            (ExtraCostMatrix *);
void                        extracostExit            (ExtraCostMatrix *);
PASTIX_INT                         extracostLoad            (ExtraCostMatrix *, FILE *);
PASTIX_INT                         extracostSave            (ExtraCostMatrix *, FILE *);

#undef static


