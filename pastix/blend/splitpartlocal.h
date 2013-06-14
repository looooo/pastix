/************************************************************/
/**                                                        **/
/**   NAME       : splitpartlocal.h                        **/
/**                                                        **/
/**   AUTHORS    : Mathieu Faverge                         **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                repartition and make processor          **/
/**                candidate groups                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 mai 2007     **/
/**                                 to     22 mai 2007     **/
/**                                                        **/
/************************************************************/

#ifndef SPLITPARTLOCAL_H
#define SPLITPARTLOCAL_H

void splitPartLocal          (BlendCtrl *, SimuCtrl *, SymbolMatrix *, const Dof *);
void propMappTreeLocal       (BlendCtrl *, const SymbolMatrix *symbmtx, const SimuCtrl *, double *);
void propMappSubtreeLocalNC  (BlendCtrl *, const SymbolMatrix *symbmtx, const SimuCtrl *, PASTIX_INT, PASTIX_INT, PASTIX_INT, double *);
void propMappSubtreeLocalOn1P(BlendCtrl *, const SymbolMatrix *symbmtx, const SimuCtrl *, PASTIX_INT, PASTIX_INT);

#endif
