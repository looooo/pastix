#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "common.h"
#include "dof.h"
#include "extendVector.h"
#include "elimin.h"
#include "cost.h"
/* #include "ftgt.h" */
#include "symbol.h"
/* #include "updown.h" */
/* #include "solver.h" */
#include "cand.h"
#include "queue.h"
#include "bulles.h"
/* #include "simu.h" */
#include "param_blend.h"
#include "blendctrl.h"
#include "assembly.h"
#include "assemblyGener.h"

void assemblyGener(pastix_int_t clustnum, Assembly1D *assemb1D, Assembly2D *assemb2D,
                   pastix_int_t clustnbr, const SymbolMatrix *symbmtx, const pastix_int_t *blprtab,
                   BlendCtrl *ctrl, const Dof * const dofptr)
/**************************************************************************/
/*  Function to gener assembly structure                                  */
/* IN SMP procnbr == clustnbr ?? (modifier blprtab dans ce fichier)       */
/* ATTENTION blprtab --> blclusttab                                       */
/**************************************************************************/
{
    pastix_int_t i, j, p;
    pastix_int_t *localnbr       = NULL;
    pastix_int_t *proc2rownbr    = NULL;
    pastix_int_t pr;
    pastix_int_t maxrow;
    double nlcoefnbr; /** Non local coeff nbr in 2D distribution / assem 1D distribution**/
    double nlbloknbr; /** Non local block nbr in 2D distribution / assem 1D distribution**/
    double lcoefnbr; /** Local coeff nbr in 2D distribution / assem 1D distribution**/
    double lbloknbr; /** Local block nbr in 2D distribution / assem 1D distribution**/
    pastix_int_t *cbprtab        = NULL;    /** cblk (in global num) to processor owner in 1D distribution**/

    pastix_int_t *bloklocalnum1D = NULL; /** Global to local 1D blocknum **/
    pastix_int_t *cblklocalnum1D = NULL; /** Global to local 2D cblknum if cblk is local in the 2D distribution **/
    pastix_int_t *cblklocalnum2D = NULL; /** Global to local 2D cblknum if cblk is local in the 2D distribution **/
    pastix_int_t *bloklocalnum2D = NULL; /** Global to local 2D blocknum **/
    pastix_int_t clustid;
    pastix_int_t bloknbr1D;
    pastix_int_t cblknbr1D;
    pastix_int_t bloknbr2D;
    pastix_int_t *bloknum        = NULL;
    pastix_int_t *cblknum        = NULL;
    pastix_int_t *blcltab        = NULL; /** blcltab[i] == Cluster owner for block i **/
    pastix_int_t delta;


    MALLOC_INTERN(localnbr, clustnbr, pastix_int_t);
    bzero(localnbr, clustnbr*sizeof(pastix_int_t));

    /** Convert blprtab in blcltab **/
    MALLOC_INTERN(blcltab, symbmtx->bloknbr, pastix_int_t);
    for(i=0;i<symbmtx->bloknbr;i++)
      {
#ifdef DEBUG_BLEND
        ASSERT(blprtab[i]>=0,MOD_BLEND);
        ASSERT(blprtab[i]<ctrl->procnbr,MOD_BLEND);
#endif
        blcltab[i] = ctrl->proc2clust[blprtab[i]];
        /*fprintf(stderr, "%ld : blprtab %ld  blcltab %ld \n", (long)i, (long)blprtab[i], (long)blcltab[i]);*/
      }

    /** Fill cbprtab **/
    MALLOC_INTERN(cbprtab, symbmtx->cblknbr, pastix_int_t);

    /***************************************************************************************************************/
    /*  if a cblk is mapped with a 1D distribution the processor owner is the same than the one decided in the     */
    /*  factorization distribution                                                                                 */
    /*  if a cblk is mapped with a 2D distribution then the processor owner in assembly is the processor that      */
    /* owns the biggest area of block in the cblk                                                                  */
    /***************************************************************************************************************/
    MALLOC_INTERN(proc2rownbr, clustnbr, pastix_int_t);
    for(i=0;i<symbmtx->cblknbr;i++)
      {
        if(ctrl->candtab[i].distrib == D1)
          /*cbprtab[i] = ctrl->proc2clust[blprtab[symbmtx->cblktab[i].bloknum]];*/
          cbprtab[i] = blcltab[symbmtx->cblktab[i].bloknum];
        else
          {
            bzero(proc2rownbr, sizeof(pastix_int_t)*clustnbr);
            for(j=symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
              /*proc2rownbr[blprtab[j]] += symbmtx->bloktab[j].lrownum - symbmtx->bloktab[j].frownum + 1;*/
              proc2rownbr[blcltab[j]] += symbmtx->bloktab[j].lrownum - symbmtx->bloktab[j].frownum + 1;

            /** Find the processor that has the largest number of row in this cblk **/
            maxrow = -1;
            pr = -1;
            for(p=0;p<clustnbr;p++)
              if(proc2rownbr[p] > maxrow)
                {
                  maxrow = proc2rownbr[p];
                  pr = p;
                }
            cbprtab[i] = pr;
          }
      }
    memFree(proc2rownbr);


    /** Compute the local numbering for cblk, blok, for each proc **/
#ifdef DEBUG_M
    ASSERT(clustnbr > 0,MOD_BLEND);
#endif
    MALLOC_INTERN(bloknum, clustnbr, pastix_int_t);
    MALLOC_INTERN(cblknum, clustnbr, pastix_int_t);

    /** GLOBAL TO LOCAL ORDERING 1D **/
    MALLOC_INTERN(cblklocalnum1D, symbmtx->cblknbr, pastix_int_t);
    MALLOC_INTERN(bloklocalnum1D, symbmtx->bloknbr, pastix_int_t);

    bzero(bloknum, sizeof(pastix_int_t)*clustnbr);
    bzero(cblknum, sizeof(pastix_int_t)*clustnbr);
    for(i=0;i<symbmtx->cblknbr;i++)
      {
        clustid = cbprtab[i];
        cblklocalnum1D[i] = cblknum[clustid];
        cblknum[clustid]++;
        for(j = symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
          {
            bloklocalnum1D[j] = bloknum[clustid];
            bloknum[clustid]++;
          }
      }
    bloknbr1D = bloknum[clustnum];
    cblknbr1D = cblknum[clustnum];

    /** GLOBAL TO LOCAL ORDERING 2D **/
    MALLOC_INTERN(cblklocalnum2D, symbmtx->cblknbr, pastix_int_t);
    MALLOC_INTERN(bloklocalnum2D, symbmtx->bloknbr, pastix_int_t);


    bzero(bloknum, sizeof(pastix_int_t)*clustnbr);
    bzero(cblknum, sizeof(pastix_int_t)*clustnbr);
    for(i=0;i<symbmtx->cblknbr;i++)
      {
        pastix_int_t flag;
        flag = 0;

        for(j = symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
          {
            /*clustid = blprtab[j];*/
            clustid = blcltab[j];
            if(clustid == clustnum)
              flag = 1; /** This is a local cblk in the 2D distribution on this processor **/

            bloklocalnum2D[j] = bloknum[clustid];
            bloknum[clustid]++;
          }
        if(flag == 1)
          {
            cblklocalnum2D[i] = cblknum[clustnum];
            cblknum[clustnum]++;
          }
        else
          cblklocalnum2D[i] = -1;
      }
    bloknbr2D = bloknum[clustnum];

    memFree(cblknum);
    memFree(bloknum);


    /** Fill assemb1D->blprtab:  ATTENTION for the 1D distribution blprtab means cbprtab **/
    MALLOC_INTERN(assemb1D->blprtab, symbmtx->cblknbr, pastix_int_t);
    memcpy(assemb1D->blprtab, cbprtab, sizeof(pastix_int_t)*symbmtx->cblknbr);


    /** Fill nocbtab **/
    MALLOC_INTERN(assemb1D->nocbtab, symbmtx->nodenbr, pastix_int_t);
    for(i=0;i<symbmtx->cblknbr;i++)
        for(j=symbmtx->cblktab[i].fcolnum;j<=symbmtx->cblktab[i].lcolnum;j++)
            assemb1D->nocbtab[j] = i;
#ifdef DEBUG_BLEND
    ASSERT(assemb1D->nocbtab[0] == 0,MOD_BLEND);
    ASSERT(assemb1D->nocbtab[symbmtx->nodenbr-1] == symbmtx->cblknbr-1,MOD_BLEND);
#endif

    /** Fill rnumtab **/
    MALLOC_INTERN(assemb1D->rnumtab, symbmtx->cblknbr, pastix_int_t);
    for(i=0;i<symbmtx->cblknbr;i++)
      {
        assemb1D->rnumtab[i] = localnbr[cbprtab[i]];
        localnbr[cbprtab[i]]++;
      }
    memFree(localnbr);

    if(ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
      {
        /*** Estimated amount of non local block ***/
        nlcoefnbr = 0;
        nlbloknbr = 0;
        lcoefnbr = 0;
        lbloknbr = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
          {
            if(ctrl->candtab[i].distrib != D1)
              {
                delta = symbmtx->cblktab[i].lcolnum - symbmtx->cblktab[i].lcolnum + 1;
                for(j=symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
                  {

                    /*if(blprtab[j] != cbprtab[i])*/
                    if(blcltab[j] != cbprtab[i])
                      {
                        nlcoefnbr += (symbmtx->bloktab[j].lrownum - symbmtx->bloktab[j].frownum + 1)*delta;
                        nlbloknbr++;
                      }
                    else
                      {
                        lcoefnbr += (symbmtx->bloktab[j].lrownum - symbmtx->bloktab[j].frownum + 1)*delta;
                        lbloknbr++;
                      }
                  }
              }

          }

        if(lbloknbr>0)
          {
            fprintf(stdout, "In assembly: in 2D distributed column block:  %g percent coef are not local \n", nlcoefnbr*100/lcoefnbr);
            fprintf(stdout, "In assembly: in 2D distributed column block:  %g percent block are not local \n", nlbloknbr*100/lbloknbr);
          }
      }


    /**************************************************************************************************/
    /*** Generation of the assembly structure 2D                                                   ****/
    /**************************************************************************************************/
    MALLOC_INTERN(assemb2D->blok2proc_tab, bloknbr1D, pastix_int_t);  /*+ local block i in 1D --> processor owner in 2D distribution +*/
    MALLOC_INTERN(assemb2D->blok2cblk_tab, bloknbr2D, pastix_int_t);  /*+ local block i in 2D --> local cblk on the same processor in
                                                                                         the 2D distribution  +*/
    MALLOC_INTERN(assemb2D->blok2blok_tab, bloknbr1D, pastix_int_t);  /*+ local block i in 1D --> local block i in the 2D distribution +*/

    for(i=0;i<symbmtx->cblknbr;i++)
      {
        if(cbprtab[i] == clustnum)
          for(j = symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
            {
              /*assemb2D->blok2proc_tab[bloklocalnum1D[j]] = blprtab[j]; */
              assemb2D->blok2proc_tab[bloklocalnum1D[j]] = blcltab[j];

              assemb2D->blok2blok_tab[bloklocalnum1D[j]] = bloklocalnum2D[j];
            }

        for(j = symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
          /*if(blprtab[j] == clustnum)*/
          if(blcltab[j] == clustnum)
            assemb2D->blok2cblk_tab[bloklocalnum2D[j]] = cblklocalnum2D[i];

      }
    /* Generation of the local symbol matrix distributed by column block for the assembly phase      */
    MALLOC_INTERN(assemb2D->symbmtx, 1, SymbolMatrix);
    symbolInit(assemb2D->symbmtx);
    symbolGener(clustnum, cblklocalnum1D, bloknbr1D, cblknbr1D, cbprtab, symbmtx, assemb2D->symbmtx, dofptr);

    memFree(cbprtab);
    memFree(blcltab);
    memFree(cblklocalnum1D);
    memFree(cblklocalnum2D);
    memFree(bloklocalnum1D);
    memFree(bloklocalnum2D);
}

void symbolGener(pastix_int_t clustnum,
                 const pastix_int_t *cblklocalnum1D,
                 pastix_int_t bloknbr1D,
                 pastix_int_t cblknbr1D,
                 const pastix_int_t *cbprtab,
                 const SymbolMatrix *symbmtx,
                 SymbolMatrix *symb1D,
                 const Dof * const dofptr)
{
  pastix_int_t i, j;
  pastix_int_t cblknum, bloknum;
  pastix_int_t nodenbr;

  symb1D->cblknbr = cblknbr1D;
  symb1D->bloknbr = bloknbr1D;

  /*************************/
  /**   Fill symb1D       **/
  /*************************/
  /* Allocations */
  MALLOC_INTERN(symb1D->cblktab, symb1D->cblknbr+1, SymbolCblk);
  MALLOC_INTERN(symb1D->bloktab, symb1D->bloknbr, SymbolBlok);

  cblknum = 0;
  bloknum = 0;
  nodenbr = 0;
  for(i=0;i<symbmtx->cblknbr;i++)
    {
      if(cbprtab[i] == clustnum)
        {
          symb1D->cblktab[cblknum].fcolnum = symbmtx->cblktab[i].fcolnum * dofptr->noddval;
          symb1D->cblktab[cblknum].lcolnum = symbmtx->cblktab[i].lcolnum * dofptr->noddval + dofptr->noddval-1;
          symb1D->cblktab[cblknum].bloknum = bloknum;
          nodenbr += symbmtx->cblktab[i].lcolnum - symbmtx->cblktab[i].fcolnum + 1;
          cblknum++;

          for(j=symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
            {
              symb1D->bloktab[bloknum].frownum = symbmtx->bloktab[j].frownum * dofptr->noddval;
              symb1D->bloktab[bloknum].lrownum = symbmtx->bloktab[j].lrownum * dofptr->noddval + dofptr->noddval-1;
              symb1D->bloktab[bloknum].cblknum = cblklocalnum1D[symbmtx->bloktab[j].cblknum];
              bloknum ++;
            }
        }
    }

  symb1D->nodenbr = nodenbr;
#ifdef DEBUG_BLEND
  ASSERT(symb1D->cblknbr == cblknum,MOD_BLEND);
  if(symb1D->bloknbr != bloknum)
    fprintf(stderr, "bloknbr %ld bloknum %ld \n", (long)symb1D->bloknbr, (long)bloknum);
  ASSERT(symb1D->bloknbr == bloknum,MOD_BLEND);
#endif

  /*  virtual cblk to avoid side effect in the loops on cblk bloks */
  symb1D->cblktab[cblknum].fcolnum = symb1D->cblktab[cblknum-1].lcolnum+1;
  symb1D->cblktab[cblknum].lcolnum = symb1D->cblktab[cblknum-1].lcolnum+1;
  symb1D->cblktab[cblknum].bloknum = bloknum;

}





#ifdef OLD_ASSEMBLY
/** AssemblyGener have to be used before solverMatrxGen (no expansion in ddl) **/
void assemblyGener(Assembly1D *assemb, pastix_int_t procnbr, const SymbolMatrix *symbmtx, const pastix_int_t *cbprtab)
/**************************************************************************/
/*  Function that was used in 1D to gener assmebly structure              */
/**************************************************************************/
{
    pastix_int_t i, j;
    pastix_int_t *localnbr;

    MALLOC_INTERN(localnbr, procnbr, pastix_int_t);
    bzero(localnbr, procnbr*sizeof(pastix_int_t));


    /** Fill blprtab **/
    /** IN 1D BL MEANS CBL **/
    MALLOC_INTERN(assemb->blprtab, symbmtx->cblknbr, pastix_int_t);
    for(i=0;i<symbmtx->cblknbr;i++)
      assemb->blprtab[i] = cbprtab[i];

    /** Fill nocbtab **/
    MALLOC_INTERN(assemb->nocbtab, symbmtx->nodenbr, pastix_int_t);
    for(i=0;i<symbmtx->cblknbr;i++)
        for(j=symbmtx->cblktab[i].fcolnum;j<=symbmtx->cblktab[i].lcolnum;j++)
            assemb->nocbtab[j] = i;

    /** Fill rnumtab **/
    MALLOC_INTERN(assemb->rnumtab, symbmtx->cblknbr, pastix_int_t);
    for(i=0;i<symbmtx->cblknbr;i++)
      {
        assemb->rnumtab[i] = localnbr[cbprtab[i]];
        localnbr[cbprtab[i]]++;
      }
    memFree(localnbr);

}
#endif
