/**
 * ftgt.h -- These lines are the data declarations
 *              for the graph structure.
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
/************************************************************/
/**                                                        **/
/**   NAME       : ftgt.h                                  **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the graph structure.                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 02 may 1998     **/
/**                                 to     16 oct 1998     **/
/**                                                        **/
/************************************************************/

#ifndef FTGT_H
#define FTGT_H

/*
 **  The type and structure definitions.
 */

/*
 ** WARNING : All structures should have a odd number of integer for memory alignment
 */

#ifndef FTGT_MAXINFO
#define FTGT_MAXINFO
/*+ Fanintarget info type +*/
typedef enum {
    FTGT_CTRBNBR = 0,           /*+ Number of contributions           +*/
    FTGT_CTRBCNT,               /*+ Number of contributions remaining +*/
    FTGT_PROCDST,               /*+ Destination for fanintarget       +*/
    FTGT_TASKDST,               /*+ Task  destination                 +*/
    FTGT_BLOKDST,               /*+ Block destination (->COMP_1D)     +*/
    FTGT_PRIONUM,               /*+ Fanintarget priority              +*/
    FTGT_FCOLNUM,               /*+ Fanintarget first column          +*/
    FTGT_LCOLNUM,               /*+ Fanintarget last column           +*/
    FTGT_FROWNUM,               /*+ Fanintarget first row             +*/
    FTGT_LROWNUM,               /*+ Fanintarget last row              +*/
#if (defined OOC) || (defined TRACE_SOPALIN) || (defined PASTIX_WITH_STARPU)
    FTGT_GCBKDST,               /*+ Global Cblk destination(->COMP_1D)+*/
    FTGT_IDTRACE,               /*+ To have 12 integer in FanInTarget +*/
#endif
    MAXINFO
} FanInInfo;
#endif /* FTGT_MAXINFO */
/*+ Fanintarget structure +*/

typedef struct FanInTarget_ {
    pastix_int_t   infotab[MAXINFO]; /*+ Fanintarget descriptor (size MAXINFO) +*/
    void          *coeftab;          /*+ Fanintarget vector access             +*/
} FanInTarget;

#endif /* FTGT_H */
