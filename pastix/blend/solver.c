/**
 * solver.c -- SolverMatrix description.
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/

#include "common.h"
#include "solver.h"

/**
 *  Function: sizeofsolver
 *
 *  Computes the size in memory of the SolverMatrix.
 *
 *  Parameters:
 *    solvptr - address of the SolverMatrix
 *
 *  Returns:
 *    SolverMatrix size.
 */
pastix_int_t
sizeofsolver(const SolverMatrix *solvptr,
                   pastix_int_t *iparm )
{
  pastix_int_t result=sizeof(SolverMatrix);
  pastix_int_t iter;
  (void)iparm;

  /* cblk and blocks arrays */
  result += solvptr->cblknbr*sizeof(SolverCblk);
  result += solvptr->bloknbr*sizeof(SolverBlok);

  /* fanin target */
  result += solvptr->ftgtnbr*sizeof(FanInTarget);
  result += solvptr->indnbr *sizeof(pastix_int_t);

  /* TODO: Check that it is not bubbletree + bubblnbr * bubblenode */
  result += solvptr->bublnbr*sizeof(BubbleTree);

  /* task */
  result += solvptr->tasknbr*sizeof(Task);

  /* ttsktab */
  result += solvptr->thrdnbr*sizeof(pastix_int_t);
  result += solvptr->thrdnbr*sizeof(pastix_int_t*);
  for (iter=0; iter<solvptr->thrdnbr; iter++)
    {
      result += solvptr->ttsknbr[iter]*sizeof(pastix_int_t);
    }

  /* proc2clust */
  result += solvptr->procnbr*sizeof(pastix_int_t);

  /* UpDownVector */
  /* TODO: 2D    */
  /* if (iparm[IPARM_DISTRIBUTION_LEVEL] == 0) */
  /*   { */
  /*     /\* UpDownCblk *\/ */
  /*     result += solvptr->cblknbr      *sizeof(UpDownCblk); */
  /*     for (iter=0; iter<solvptr->cblknbr; iter++) */
  /*       { */
  /*         /\* browproctab / browcblktab *\/ */
  /*         result += 2*solvptr->updovct.cblktab[iter].browprocnbr*sizeof(pastix_int_t); */
  /*       } */
  /*     /\* gcblk2list *\/ */
  /*     result += solvptr->updovct.gcblk2listnbr*sizeof(pastix_int_t); */
  /*     /\* listptr *\/ */
  /*     result += solvptr->updovct.listptrnbr   *sizeof(pastix_int_t); */
  /*     /\* listcblk / listblok *\/ */
  /*     result += 2*solvptr->updovct.listnbr    *sizeof(pastix_int_t); */
  /*     /\* loc2glob *\/ */
  /*     result += solvptr->updovct.loc2globnbr  *sizeof(pastix_int_t); */
  /*     /\* lblk2gcblk *\/ */
  /*     result += solvptr->bloknbr              *sizeof(pastix_int_t); */
  /*   } */

  return result;
}
