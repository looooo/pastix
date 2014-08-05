/**
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
#include "z_raff_functions.h"

/*
 ** Section: Functions declarations
 */

/* Raffinement du second membre */
#define z_bicgstab_smp         API_CALL(z_bicgstab_smp)
#define z_bicgstab_thread      API_CALL(z_bicgstab_thread)

void* z_bicgstab_smp(void *arg);

/* Lancement d'une des fonctions seules */
void z_bicgstab_thread(z_SolverMatrix *datacode, z_SopalinParam *sopaparam);


/*
 ** Section: Threads routines
 */

/*
 * Function: API_CALL(z_bicgstab_smp)
 *
 * Refine the solution.
 *
 * Computes :
 *
 * Algorithm taken from Saad Iterative Methods' Book.
 * Second Edition, p219
 *
 * Parameters:
 *   arg - Pointer to a <sopthread_data_t> structure containing
 *         the <z_Sopalin_Data_t> structure and the thread number ID.
 */
void* API_CALL(z_bicgstab_smp) ( void *arg )
{
  /* Choix du solveur */
  struct z_solver solveur = {NULL};
  z_Pastix_Solveur(&solveur);

  Clock   raff_clk;
  double  t0      = 0;
  double  t3      = 0;
  z_RAFF_INT     n       = solveur.N(arg);
  z_RAFF_FLOAT   normb   = 0.0;
  z_RAFF_FLOAT   normr   = 0.0;
  int     nb_iter = 0;
  z_RAFF_FLOAT   epsilon = solveur.Eps(arg);
  z_RAFF_INT     itermax = solveur.Itermax(arg);
  z_RAFF_FLOAT   tmp     = 0.0;

  z_RAFF_FLOAT * gradb  = NULL; /* Second membre b */
  z_RAFF_FLOAT * gradr  = NULL; /* Solution actuelle */
  z_RAFF_FLOAT * gradr2 = NULL; /* Condition initiale bis r^ */
  z_RAFF_FLOAT * gradp  = NULL;
  z_RAFF_FLOAT * grady  = NULL;
  z_RAFF_FLOAT * gradv  = NULL;
  z_RAFF_FLOAT * grads  = NULL;
  z_RAFF_FLOAT * gradz  = NULL;
  z_RAFF_FLOAT * gradt  = NULL;
  z_RAFF_FLOAT * grad2  = NULL; /* Vecteurs de transition */
  z_RAFF_FLOAT * grad3  = NULL;

  /* Alpha et Beta ne sont utilisés que par le thread 0 */
  z_RAFF_FLOAT * alpha = NULL;
  z_RAFF_FLOAT * beta  = NULL;
  z_RAFF_FLOAT * v1    = NULL;
  z_RAFF_FLOAT * v2    = NULL;
  z_RAFF_FLOAT * w     = NULL;
  z_RAFF_FLOAT * gradx = NULL;

  gradb  = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));
  gradr  = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));
  gradr2 = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));
  gradp  = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));
  grady  = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));
  gradv  = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));
  grads  = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));
  gradz  = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));
  gradt  = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));
  grad2  = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));
  grad3  = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));
  alpha  = solveur.Malloc(arg, 1 * sizeof(z_RAFF_FLOAT));
  beta   = solveur.Malloc(arg, 1 * sizeof(z_RAFF_FLOAT));
  v1     = solveur.Malloc(arg, 1 * sizeof(z_RAFF_FLOAT));
  v2     = solveur.Malloc(arg, 1 * sizeof(z_RAFF_FLOAT));
  w      = solveur.Malloc(arg, 1 * sizeof(z_RAFF_FLOAT));
  gradx  = solveur.Malloc(arg, n * sizeof(z_RAFF_FLOAT));

  gradb  = solveur.Synchro(arg, (void*) gradb,  0);
  gradr  = solveur.Synchro(arg, (void*) gradr,  1);
  gradr2 = solveur.Synchro(arg, (void*) gradr2, 2);
  gradp  = solveur.Synchro(arg, (void*) gradp,  3);
  grady  = solveur.Synchro(arg, (void*) grady,  4);
  gradv  = solveur.Synchro(arg, (void*) gradv,  5);
  grads  = solveur.Synchro(arg, (void*) grads,  6);
  gradz  = solveur.Synchro(arg, (void*) gradz,  7);
  gradt  = solveur.Synchro(arg, (void*) gradt,  8);
  grad2  = solveur.Synchro(arg, (void*) grad2,  9);
  grad3  = solveur.Synchro(arg, (void*) grad3,  10);
  gradx  = solveur.Synchro(arg, (void*) gradx,  11);

  alpha  = solveur.Synchro(arg, (void*) alpha,  12);
  beta   = solveur.Synchro(arg, (void*) beta,   13);
  w      = solveur.Synchro(arg, (void*) w,      14);

  RAFF_CLOCK_INIT;

  solveur.B(arg, gradb);
  solveur.X(arg, gradx);

  /* r = b - Ax */
  solveur.bMAx(arg, gradb, gradx, gradr);
  normb = solveur.Norm(arg, gradb);
  normr = solveur.Norm(arg, gradr);

  /* r2 = r */
  solveur.Copy(arg, gradr, gradr2, 0);
  /* p = r */
  solveur.Copy(arg, gradr, gradp, 1);

  /* tmp = ||r|| / ||b|| */
  tmp = normr / normb;

  while (((float)tmp > (float)epsilon) && (nb_iter < itermax))
    {
      RAFF_CLOCK_STOP;
      t0 = RAFF_CLOCK_GET;

      nb_iter++;

      /* y = M-1 * p */
      solveur.Precond(arg, gradp, grady, 1);

      /* v = Ay */
      solveur.Ax(arg, grady, gradv);

      /* alpha = (r, r2) / (v, r2) */
      /* alpha = (v, r2) */
      solveur.Dotc(arg, gradv, gradr2, alpha, 0);
      /* beta = (r, r2) */
      solveur.Dotc(arg, gradr, gradr2, beta, 1);

      /* alpha = beta / alpha : alpha = (r, r2) / (v, r2) */
      solveur.Div(arg, beta, alpha, alpha, 0);

      /* s = r - alpha * v */
      solveur.Copy(arg, gradr, grads, 1);
      solveur.AXPY(arg, -1, alpha, grads, gradv, 1);

      /* z = M-1s */
      solveur.Precond(arg, grads, gradz, 1);

      /* t = Az */
      solveur.Ax(arg, gradz, gradt);

      /* w = (M-1t, M-1s) / (M-1t, M-1t) */
      /* grad2 = M-1t */
      solveur.Precond(arg, gradt, grad2, 1);

      /* v1 = (M-1t, M-1s) */
      /* v2 = (M-1t, M-1t) */
      solveur.Dotc(arg, gradz, grad2, v1, 0);
      solveur.Dotc(arg, grad2, grad2, v2, 1);

      solveur.Div(arg, v1, v2, w, 1);

      /* x = x + alpha * y + w * z */
      /* x = x + alpha * y */
      solveur.AXPY(arg, 1, alpha, gradx, grady, 0);

      /* x = x + w * z */
      solveur.AXPY(arg, 1, w, gradx, gradz, 0);

      /* r = s - w * t*/
      solveur.Copy(arg, grads, gradr, 1);
      solveur.AXPY(arg, -1, w, gradr, gradt, 1);

      /* beta = (r', r2) / (r, r2) * (alpha / w) */
      /* v1 = (r', r2) */
      solveur.Dotc(arg, gradr, gradr2, v1, 1);

      /* v2 = alpha / w */
      solveur.Div(arg, alpha, w, v2, 0);

      /* beta = v1 / beta */
      solveur.Div(arg, v1, beta, beta, 0);

      /* beta = beta * v2 */
      solveur.Mult(arg, beta, v2, beta, 1);

      /* p = r + beta * (p - w * v) */
      /* p = p - w * v */
      solveur.AXPY(arg, -1, w, gradp, gradv, 1);

      /* p = r + beta * p */
      solveur.BYPX(arg, beta, gradr, gradp, 1);

      normr = solveur.Norm(arg, gradr);

      RAFF_CLOCK_STOP;
      t3 = RAFF_CLOCK_GET;

      tmp = normr / normb;
      solveur.Verbose(arg, t0, t3, tmp, nb_iter);
    }

  solveur.End(arg, tmp, nb_iter, t3, gradx);

  solveur.Free(arg, (void*) gradb);
  solveur.Free(arg, (void*) gradr);
  solveur.Free(arg, (void*) gradr2);
  solveur.Free(arg, (void*) gradp);
  solveur.Free(arg, (void*) grady);
  solveur.Free(arg, (void*) gradv);
  solveur.Free(arg, (void*) grads);
  solveur.Free(arg, (void*) gradz);
  solveur.Free(arg, (void*) gradt);
  solveur.Free(arg, (void*) grad2);
  solveur.Free(arg, (void*) grad3);
  solveur.Free(arg, (void*) alpha);
  solveur.Free(arg, (void*) beta);
  solveur.Free(arg, (void*) v1);
  solveur.Free(arg, (void*) v2);
  solveur.Free(arg, (void*) w);
  solveur.Free(arg, (void*) gradx);
  return 0;
}


/*
** Section: Function creating threads
*/
void API_CALL(z_bicgstab_thread)(z_SolverMatrix *datacode, z_SopalinParam *sopaparam)
{
  z_raff_thread(datacode, sopaparam, &API_CALL(z_bicgstab_smp));
}
