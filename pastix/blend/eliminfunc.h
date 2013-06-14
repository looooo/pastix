/************************************************************/
/**                                                        **/
/**   NAME       : eliminfunc.h                            **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Functions using the elimination         **/
/**                structures                              **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     27 jul 1998     **/
/**                                                        **/
/************************************************************/
#ifndef ELIMINFUNC_H
#define static
#endif

void                        eliminGraphBuild  (const SymbolMatrix *, EliminGraph *);
void                        eliminTreeBuild   (const SymbolMatrix *, BlendCtrl *);
PASTIX_INT                         treeLeaveNbr      (const EliminTree *);
PASTIX_INT                         treeLevel         (const EliminTree *);
PASTIX_INT                         nodeTreeLevel     (PASTIX_INT, const EliminTree *);
static void                 etreeBuild        (EliminTree *, const SymbolMatrix *);

#undef static
