#include "common.h"
#include "symbol.h"
#include "cost.h"
#include "extrastruct.h"
#include "dof.h"
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
/* #include "splitpart.h" */
/* #include "assert.h" */
#include "splitfunc.h"


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
