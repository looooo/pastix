/**
 * z_solver.c -- z_SolverMatrix description.
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/

#include "common.h"
#include "z_solver.h"

/**
 *  Function: sizeofsolver
 *
 *  Computes the size in memory of the z_SolverMatrix.
 *
 *  Parameters:
 *    solvptr - address of the z_SolverMatrix
 *
 *  Returns:
 *    z_SolverMatrix size.
 */
pastix_int_t
z_sizeofsolver(const z_SolverMatrix *solvptr,
                      pastix_int_t *iparm )
{
  pastix_int_t result=sizeof(z_SolverMatrix);
  pastix_int_t iter;

  /* cblk and blocks arrays */
  result += solvptr->cblknbr*sizeof(z_SolverCblk);
  result += solvptr->bloknbr*sizeof(z_SolverBlok);

  /* fanin target */
  result += solvptr->ftgtnbr*sizeof(z_FanInTarget);
  result += solvptr->indnbr *sizeof(pastix_int_t);

  /* TODO: Check that it is not bubbletree + bubblnbr * bubblenode */
  result += solvptr->bublnbr*sizeof(BubbleTree);

  /* task */
  result += solvptr->tasknbr*sizeof(z_Task);

  /* ttsktab */
  result += solvptr->thrdnbr*sizeof(pastix_int_t);
  result += solvptr->thrdnbr*sizeof(pastix_int_t*);
  for (iter=0; iter<solvptr->thrdnbr; iter++)
    {
      result += solvptr->ttsknbr[iter]*sizeof(pastix_int_t);
    }

  /* proc2clust */
  result += solvptr->procnbr*sizeof(pastix_int_t);

  /* z_UpDownVector */
  /* TODO: 2D    */
  if (iparm[IPARM_DISTRIBUTION_LEVEL] == 0)
    {
      /* UpDownCblk */
      result += solvptr->cblknbr      *sizeof(z_UpDownCblk);
      for (iter=0; iter<solvptr->cblknbr; iter++)
        {
          /* browproctab / browcblktab */
          result += 2*solvptr->updovct.cblktab[iter].browprocnbr*sizeof(pastix_int_t);
        }
      /* gcblk2list */
      result += solvptr->updovct.gcblk2listnbr*sizeof(pastix_int_t);
      /* listptr */
      result += solvptr->updovct.listptrnbr   *sizeof(pastix_int_t);
      /* listcblk / listblok */
      result += 2*solvptr->updovct.listnbr    *sizeof(pastix_int_t);
      /* loc2glob */
      result += solvptr->updovct.loc2globnbr  *sizeof(pastix_int_t);
      /* lblk2gcblk */
      result += solvptr->bloknbr              *sizeof(pastix_int_t);
    }

  return result;
}
