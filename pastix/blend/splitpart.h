/************************************************************/
/**                                                        **/
/**   NAME       : splitpart.h                             **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                repartition and make processor          **/
/**                candidate groups                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     09 sep 1998     **/
/**                                                        **/
/************************************************************/
#ifndef SPLITPART_H

#define static

#endif

void propMappTree( BlendCtrl         *ctrl,
                   EliminTree        *etree,
                   Cand              *candtab,
                   CostMatrix        *costmtx,
                   SymbolMatrix      *symbmtx,
                   const Dof         *dofptr,
                   ExtraSymbolMatrix *extrasymb,
                   ExtraCostMatrix   *extracost,
                   pastix_int_t       procnbr,
                   int split, int nocrossproc, int allcand );

void  splitPart     (SymbolMatrix *, BlendCtrl *, const Dof *);
pastix_int_t   check_candidat(SymbolMatrix *, BlendCtrl *);

static void setDistribType       (const pastix_int_t, SymbolMatrix *, Cand *, const pastix_int_t);
static void setSubtreeDistribType(const SymbolMatrix *, const CostMatrix *, pastix_int_t , const BlendCtrl *, pastix_int_t);

static void splitOnProcs    (SymbolMatrix *, ExtraSymbolMatrix *, ExtraCostMatrix *, BlendCtrl *, 
			     const Dof *, pastix_int_t, pastix_int_t);
static void splitCblk       (SymbolMatrix *, ExtraSymbolMatrix *, ExtraCostMatrix *, BlendCtrl *, 
			     const Dof *, pastix_int_t, pastix_int_t, pastix_int_t *);

static void printTree          (FILE*, const EliminTree *, pastix_int_t);

static double maxProcCost     (double *, pastix_int_t);
static double blokUpdateCost  (pastix_int_t, pastix_int_t, CostMatrix *, ExtraCostMatrix *, const SymbolMatrix *, 
			       const ExtraSymbolMatrix *, BlendCtrl *, const Dof *);

static pastix_int_t    countBlok            (pastix_int_t, SymbolMatrix *, pastix_int_t);
static pastix_int_t    setSubtreeBlokNbr    (pastix_int_t, const EliminTree *, SymbolMatrix *, ExtraSymbolMatrix *, pastix_int_t);
static void   clusterCandCorrect   (pastix_int_t, Cand *, const EliminTree *, BlendCtrl *);
static void   setClusterCand       (pastix_int_t, Cand *, const EliminTree *, pastix_int_t, pastix_int_t);

#undef static

