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
void       assemblyGener(Assembly1D *, pastix_int_t, const SymbolMatrix *, const pastix_int_t *);
#else
void assemblyGener(pastix_int_t clustnum, Assembly1D *assemb1D, Assembly2D *assemb2D, pastix_int_t clustnbr, const SymbolMatrix *symbmtx, const pastix_int_t *blprtab, BlendCtrl *ctrl, const Dof * const dofptr);
static void symbolGener(pastix_int_t clustnum, const pastix_int_t *cblklocalnum1D, pastix_int_t bloknbr1D, pastix_int_t cblknbr1D, const pastix_int_t *cbprtab, const SymbolMatrix *symbmtx, SymbolMatrix *symb1D, const Dof * const dofptr);
#endif

#undef static
