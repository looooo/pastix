/**
 *
 * @file extraclk.h
 *
 *  PaStiX analyse routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _EXTRACBLK_H_
#define _EXTRACBLK_H_

typedef struct ExtraCblk_s {
    pastix_int_t  cblknbr;
    pastix_int_t  addcblk;  /*+ Number of cblk created                            +*/
    pastix_int_t  addblok;  /*+ Number of blok created                            +*/
    pastix_int_t  addblof;  /*+ Number of blok created due to facing cblk splited +*/
    pastix_int_t *sptcblk;  /*+ Index for splitted cblk in the cblktab            +*/
    pastix_int_t *sptcbnb;  /*+ Number of splitted cblk for a cblk                +*/
    pastix_int_t  curcblk;  /*+ Cursor for cblktab                                +*/
    pastix_int_t  sizcblk;  /*+ Size of allocated cblktab                         +*/
    SymbolCblk   *cblktab;  /*+ Array of column blocks [+1,based]                 +*/
} ExtraCblk_t;

void extraCblkInit( pastix_int_t cblknbr,
                    ExtraCblk_t *extracblk );

void extraCblkExit( ExtraCblk_t *extracblk );

void extraCblkAdd( ExtraCblk_t *extracblk,
                   pastix_int_t fcolnum,
                   pastix_int_t lcolnum );

void extraCblkMerge( ExtraCblk_t  *extracblk,
                     SymbolMatrix *newsymb,
                     Cand        **candtab );

#endif /* _EXTRACBLK_H_ */
