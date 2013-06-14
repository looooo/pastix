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
  pastix_int_t                       baseval;              /*+ Base value for numberings                    +*/
  pastix_int_t                       addcblk;              /*+ Number of cblk created                       +*/
  pastix_int_t                       addblok;              /*+ Number of blok created                       +*/
  pastix_int_t        *              sptcblk;              /*+ Index for splitted cblk in the cblktab       +*/
  pastix_int_t        *              sptcbnb;              /*+ Number of splitted cblk for a cblk           +*/
  pastix_int_t        *              sptblok;              /*+ Index for splitted blok in the bloktab       +*/
  pastix_int_t        *              sptblnb;              /*+ Number of splitted blok for a blok           +*/
  pastix_int_t        *              subtreeblnbr;         /*+ Number of blok in the subtree                +*/
  pastix_int_t                       curcblk;              /*+ Cursor for cblktab                           +*/
  pastix_int_t                       sizcblk;              /*+ Size of allocated cblktab                    +*/
  pastix_int_t                       curblok;              /*+ Cursor for bloktab                           +*/
  pastix_int_t                       sizblok;              /*+ Size of allocated bloktab                    +*/
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
pastix_int_t                         extrasymbolInit          (ExtraSymbolMatrix *);
void                        extrasymbolExit          (ExtraSymbolMatrix *);
pastix_int_t                         extrasymbolLoad          (ExtraSymbolMatrix *, FILE *);
pastix_int_t                         extrasymbolSave          (ExtraSymbolMatrix *, FILE *);
pastix_int_t                         extracostInit            (ExtraCostMatrix *);
void                        extracostExit            (ExtraCostMatrix *);
pastix_int_t                         extracostLoad            (ExtraCostMatrix *, FILE *);
pastix_int_t                         extracostSave            (ExtraCostMatrix *, FILE *);

#undef static


