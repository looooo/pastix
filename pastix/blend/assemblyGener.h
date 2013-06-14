/************************************************************/
/**                                                        **/
/**   NAME       : assemblyGener.h                         **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Generation of the assembly strucuture   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     09 sep 1998     **/
/**                                                        **/
/************************************************************/
#ifndef ASSEMBLY
#define static
#endif

#ifdef OLD_ASSEMBLY
void       assemblyGener(Assembly1D *, PASTIX_INT, const SymbolMatrix *, const PASTIX_INT *);
#else
void assemblyGener(PASTIX_INT clustnum, Assembly1D *assemb1D, Assembly2D *assemb2D, PASTIX_INT clustnbr, const SymbolMatrix *symbmtx, const PASTIX_INT *blprtab, BlendCtrl *ctrl, const Dof * const dofptr);
static void symbolGener(PASTIX_INT clustnum, const PASTIX_INT *cblklocalnum1D, PASTIX_INT bloknbr1D, PASTIX_INT cblknbr1D, const PASTIX_INT *cbprtab, const SymbolMatrix *symbmtx, SymbolMatrix *symb1D, const Dof * const dofptr);
#endif

#undef static
