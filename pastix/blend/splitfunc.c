#include "common.h"
#include "symbol.h"
#include "cost.h"
#include "extrastruct.h"
#include "dof.h"
#include "param_blend.h"
#include "elimin.h"
#include "cand.h"
#include "queue.h"
#include "extendVector.h"
#include "bulles.h"
#include "blendctrl.h"
#include "ftgt.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "simu.h"
#include "costfunc.h"
/* #include "assembly.h" */
/* #include "eliminfunc.h" */
/* #include "splitpart.h" */
/* #include "assert.h" */
#include "splitfunc.h"

pastix_int_t splitSeqCblk2D( pastix_int_t cblknum, pastix_int_t procnbr, const SymbolMatrix * symbptr, const ExtraSymbolMatrix * extrasymbptr, const Dof * dofptr, const BlendCtrl * ctrl,
                    pastix_int_t (*P)(pastix_int_t, pastix_int_t, const SymbolMatrix *,const ExtraSymbolMatrix *, const Dof *, pastix_int_t, const BlendCtrl *) )
{
  pastix_int_t a, c;
  a = 1;
  c = symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum + 1;

  /*(!P(cblknum, procnbr, symbptr, extrasymbptr, dofptr, 1, ctrl))
    return 1;*/
  if( P(cblknum, procnbr, symbptr, extrasymbptr, dofptr, c, ctrl))
    return c;

  while(a<c-1)
    {
      if(P(cblknum, procnbr, symbptr, extrasymbptr, dofptr, (a+c)/2, ctrl))
        a = (a+c)/2;
      else
        c = (a+c)/2;
    }
#ifdef DEBUG_BLEND
  ASSERT(!P(cblknum, procnbr, symbptr, extrasymbptr, dofptr, c, ctrl),MOD_BLEND);
  ASSERT(a>=1 && a<= symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum + 1,MOD_BLEND);
#endif

  return a;
}

#if 0
/** Basic criterion based on broadness of cblk **/
pastix_int_t P1D(pastix_int_t cblknum, pastix_int_t procnbr, const SymbolMatrix *symbptr, const ExtraSymbolMatrix * const extrasymbptr, const Dof * dofptr, pastix_int_t nseq, const BlendCtrl * ctrl)
{
  pastix_int_t delta = symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum + 1;
  if( dofptr->noddval*delta/nseq >= ctrl->option->blcolmin )
    return 1;
  return 0;
}

/** Criterion based on overhead du to the splitting **/
pastix_int_t P2D(pastix_int_t cblknum, pastix_int_t procnbr, const SymbolMatrix * symbptr, const ExtraSymbolMatrix * const extrasymbptr, const Dof * dofptr, pastix_int_t nseq, const BlendCtrl * ctrl)
{
  pastix_int_t bloknbr;
  pastix_int_t sbloknbr;
  SymbolBlok *bloktab  = NULL;
  SymbolBlok *sbloktab = NULL;
  pastix_int_t        *indtab  = NULL;
  double cost1;
  double cost2;


  pastix_int_t k;

  /* Number of blocks */
  bloknbr = cblkNbr(cblknum, symbptr, extrasymbptr);

  MALLOC_INTERN(bloktab, bloknbr, SymbolBlok);
  build_cblk(cblknum, symbptr, extrasymbptr, bloktab);

  /* Number of splitted blocks */
  sbloknbr = (nseq*(nseq+1))/2 + (bloknbr-1)*nseq;


  /* Cost of the non splitted cblk */
  cost1 = ctrl->costmtx->cblktab[cblknum].total;
#ifdef DEBUG_BLEND
  if(cost1 != cblkCost(bloknbr, bloktab, dofptr))
    fprintf(stderr, "cblknum %ld ; nseq %ld bloknbr %ld costmtx %g mycost %g \n", (long)cblknum,  (long)nseq, (long)bloknbr, cost1, cblkCost(bloknbr, bloktab, dofptr));
/*  ASSERT(cost1 == cblkCost(bloknbr, bloktab, dofptr),MOD_BLEND);*/
#endif

  /** give the structure of splitted cblk (it is not repercuted on the symbol matrix) **/
  MALLOC_INTERN(indtab,   nseq+1,   pastix_int_t);
  MALLOC_INTERN(sbloktab, sbloknbr, SymbolBlok);

  virtualSplit(nseq, bloknbr, bloktab, indtab, sbloktab);

  /* Cost of splitted cblk */
  cost2 = 0;
  for(k=0;k<nseq;k++)
    {
      cost2 += cblkCost(bloknbr+nseq-k-1, &(sbloktab[indtab[k]]), dofptr);
    }

  memFree(bloktab);
  memFree(indtab);
  memFree(sbloktab);
  fprintf(stdout, "cblknum %ld NSEQ %ld COST2 : %g COST1 : %g Cost2/Cost1 %g\n",  (long)cblknum,(long)nseq, cost2, cost1, cost2/cost1);
  if(cost1==0)
    return 0;

  procnbr = MIN(procnbr, nseq/4);

  if(cost2/cost1 < 0.5*procnbr)
    return 1;
  return 0;

}
#endif

/** OIMBE TO DO: le calcul total du cout splitter peut etre bp plus rapide en groupant
  les calculs sur les odb dans les cblk decouper **/
void virtualSplit(pastix_int_t nseq, pastix_int_t bloknbr, const SymbolBlok * src_bloktab,
                         pastix_int_t *indtab, SymbolBlok * dest_bloktab)
{
  pastix_int_t k;
  pastix_int_t i;
  pastix_int_t step;

  /** Index of splitted cblk in dest_bloktab **/
  /* + a virtual extra end cblk */

  indtab[0] = 0;
  for(k=1;k<=nseq;k++)
    indtab[k] = indtab[k-1] + bloknbr + nseq-k;

  step = (src_bloktab[0].lrownum - src_bloktab[0].frownum + 1)/nseq;

  /** Compute the splitted cblk **/
  for(k=0;k<nseq;k++)
    {
      /* We just have to compute diagonal splitting */
      for(i=0;i < nseq-k-1;i++)
        {
          dest_bloktab[i+indtab[k]].frownum = src_bloktab[0].frownum + (i+k)*step;
          dest_bloktab[i+indtab[k]].lrownum = src_bloktab[0].frownum + (i+k+1)*step-1;
        }
      dest_bloktab[i+indtab[k]].frownum = src_bloktab[0].frownum + (i+k)*step;
      dest_bloktab[i+indtab[k]].lrownum = src_bloktab[0].lrownum;
      i++;
      /* The odbs remain the same */
      if(bloknbr>1)
        memcpy(&(dest_bloktab[i+indtab[k]]), &(src_bloktab[1]), sizeof(SymbolBlok)*(bloknbr-1));
    }
  /*fprintf(stdout, "CBlok Orig \n");
  for(i=0;i<bloknbr;i++)
    printf("O [%ld %ld]\n",src_bloktab[i].frownum,src_bloktab[i].lrownum );

  for(k=0;k<nseq;k++)
    {
      printf("CBlok %ld \n", k);
      for(i=indtab[k];i<indtab[k+1];i++)
        printf("S [%ld %ld]\n",dest_bloktab[i].frownum,dest_bloktab[i].lrownum );
    }*/

}

/** Assure that cblkComputeCost and cblkCost compute the same things !!!! **/
double cblkCost(pastix_int_t bloknbr, const SymbolBlok * bloktab, const Dof * dofptr)
{
    pastix_int_t l, h, g;
    pastix_int_t k;
    double total_cost = 0;
    double compute_cost = 0;
    double send_cost    = 0;
    double contrib_cost = 0;
    /** we need the height of cblk non empty lines  and the broadness
      of the cblk to compute the local compute cost **/
#ifdef DOF_CONSTANT
    l = (bloktab[0].lrownum -bloktab[0].frownum + 1)*(dofptr)->noddval;
#endif

    g = 0;
    for(k=0;k<bloknbr;k++)
      {
#ifdef  DOF_CONSTANT
        g += (bloktab[k].lrownum - bloktab[k].frownum + 1)*(dofptr)->noddval;
#endif
      }

    /** retrieve diag height so let g be the odb non empty lines height **/
    g -= l;

    /** compute the local compute cost **/
    if(l!=0)
        compute_cost += computeCost(l, g);
    else
      compute_cost = 0;

    /** compute for each odb its contribution compute cost and add cost **/
    for(k=1;k<bloknbr;k++)
      {
#ifdef  DOF_CONSTANT
        h = (bloktab[k].lrownum - bloktab[k].frownum + 1)*(dofptr)->noddval;
#endif
        /* g is the odb lines number above this odb (odb lines include)*/
        /*if(l!=0 && h != 0 && g != 0)*/
        contrib_cost     = contribCompCost(l, h, g);
        /*if(h != 0 && g != 0)*/
        contrib_cost     += contribAddCost(h, g);
        send_cost += contrib_cost;
        g -= h;
      }
    total_cost = compute_cost + send_cost;
    return total_cost;
}

pastix_int_t cblkNbr(pastix_int_t cblknum,  const SymbolMatrix * symbptr, const ExtraSymbolMatrix * extrasymbptr)
{
  pastix_int_t bloknbr = 0;
  pastix_int_t i;
  for(i=symbptr->cblktab[cblknum].bloknum;i<symbptr->cblktab[cblknum+1].bloknum;i++)
    {
      if(extrasymbptr->sptblnb[i]>1)
        bloknbr += extrasymbptr->sptblnb[i];
      else
        bloknbr++;
    }
  return bloknbr;
}

void build_cblk(pastix_int_t cblknum, const SymbolMatrix * symbptr, const ExtraSymbolMatrix * extrasymbptr, SymbolBlok *bloktab)
{
  pastix_int_t i;
  pastix_int_t blokcur = 0;
  pastix_int_t sptnbr;
  for(i=symbptr->cblktab[cblknum].bloknum;i<symbptr->cblktab[cblknum+1].bloknum;i++)
    {
      if(extrasymbptr->sptblok[i]>=0)
        {
          sptnbr = extrasymbptr->sptblnb[i];
#ifdef DEBUG_BLEND
          ASSERT(sptnbr>0,MOD_BLEND);
#endif
          memcpy(&(bloktab[blokcur]), &(extrasymbptr->bloktab[extrasymbptr->sptblok[i]]), sptnbr*sizeof(SymbolBlok));
          blokcur += sptnbr;
        }
      else
        {
          bloktab[blokcur].frownum = symbptr->bloktab[i].frownum;
          bloktab[blokcur].lrownum = symbptr->bloktab[i].lrownum;
          blokcur++;
        }
    }

}
