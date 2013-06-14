/************************************************************/
/**                                                        **/
/**   NAME       : assembly.h                              **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the assembly process.               **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 jul 1998     **/
/**                                 to     08 sep 1998     **/
/**                                                        **/
/************************************************************/

#define ASSEMBLY_H

#ifdef CXREF_DOC
#include "common.h"
#endif /* CXREF_DOC */

/*
**  The type and structure definitions.
*/

/*+ 1D Assembly structure based on the 1D distribution for assembly. +*/

typedef struct Assembly1D_ {
  pastix_int_t *                     blprtab;              /*+ Array of block to processor mapping        +*/
  pastix_int_t *                     nocbtab;              /*+ Array of node to column block (1D) mapping +*/
  pastix_int_t *                     rnumtab;              /*+ Absolute number to relative number tab     +*/
} Assembly1D;


/*+ 2D Assembly structure. +*/
typedef struct Assembly2D_ {
  SymbolMatrix *            symbmtx;              /*+ Symbol matrix in 1D distribution needed for assembling the matrix   +*/ 
  pastix_int_t *                     blok2proc_tab;        /*+ local block i in 1D --> processor owner in 2D distribution +*/
  pastix_int_t *                     blok2cblk_tab;        /*+ local block i in 2D --> local cblk on the same processor in the 2D distribution  +*/
  pastix_int_t *                     blok2blok_tab;        /*+ local block i in 1D --> local block i in the 2D distribution +*/
} Assembly2D;
