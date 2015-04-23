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

void* z_grad_smp         (void *arg);

/* Lancement d'une des fonctions seules */
void z_grad_thread (z_SolverMatrix *datacode, z_SopalinParam *sopaparam);
/*
 ** Section: Threads routines
 */


/*
 Function: z_grad_smp

 Refine the solution using conjugate gradian method.

 Parameters:
 arg - Pointer to a <sopthread_data_t> structure containing
 the <z_Sopalin_Data_t> structure and the thread number ID.

 */
void* z_grad_smp(void *arg)
{
  /* Choix du solveur */
  struct solver solveur = {NULL};
  Pastix_Solveur(&solveur);

  /* Variables */
  Clock  raff_clk;
  double t0      = 0;
  double t3      = 0;
  RAFF_FLOAT  tmp     = 0.0;
  RAFF_FLOAT  normr;
  RAFF_FLOAT  normb;
  RAFF_FLOAT  epsilon = solveur.Eps(arg);
  RAFF_INT    itermax = solveur.Itermax(arg);
  RAFF_INT    nb_iter = 0;
  RAFF_INT    n       = solveur.N(arg);

  RAFF_FLOAT *gradb = NULL;
  RAFF_FLOAT *gradr = NULL;
  RAFF_FLOAT *gradp = NULL;
  RAFF_FLOAT *gradz = NULL;
  RAFF_FLOAT *grad2 = NULL;
  RAFF_FLOAT *alpha = NULL;
  RAFF_FLOAT *beta  = NULL;
  RAFF_FLOAT *gradx = NULL;

  /* Initialisation des vecteurs */
  gradb = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));
  gradr = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));
  gradp = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));
  gradz = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));
  grad2 = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));
  alpha = solveur.Malloc(arg,     sizeof(RAFF_FLOAT));
  beta  = solveur.Malloc(arg,     sizeof(RAFF_FLOAT));
  gradx = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));

  gradb = solveur.Synchro(arg, (void*) gradb, 0);
  gradr = solveur.Synchro(arg, (void*) gradr, 1);
  gradp = solveur.Synchro(arg, (void*) gradp, 2);
  gradz = solveur.Synchro(arg, (void*) gradz, 3);
  grad2 = solveur.Synchro(arg, (void*) grad2, 4);
  gradx = solveur.Synchro(arg, (void*) gradx, 5);
  alpha = solveur.Synchro(arg, (void*) alpha, 6);
  beta  = solveur.Synchro(arg, (void*) beta,  7);

  RAFF_CLOCK_INIT;

  solveur.B(arg, gradb);
  solveur.X(arg, gradx);

  /* r=b-ax */
  solveur.bMAx(arg, gradb, gradx, gradr);
  normb = solveur.Norm(arg, gradb);
  normr = solveur.Norm(arg, gradr);
  tmp = normr / normb;

  /* z = M-1 r */
  solveur.Precond(arg, gradr, gradz, 1);

  solveur.Copy(arg, gradz, gradp, 1);

  while (((double)tmp > (double)epsilon) && (nb_iter < itermax))
    {
      RAFF_CLOCK_STOP;
      t0 = RAFF_CLOCK_GET;
      nb_iter++;

      /* grad2 = A * p */
      solveur.Ax(arg, gradp, grad2);

      /* alpha = <r, z> / <Ap, p> */
      solveur.Dotc(arg, gradr, gradz, beta, 0);
      solveur.Dotc(arg, grad2, gradp, alpha, 1);
      solveur.Div(arg, beta, alpha, alpha, 1);

      /* x = x + alpha * p */
      solveur.AXPY(arg, 1, alpha, gradx, gradp, 0);

      /* r = r - alpha * A * p */
      solveur.AXPY(arg, -1, alpha, gradr, grad2, 1);

      /* z = M-1 * r */
      solveur.Precond(arg, gradr, gradz, 1);

      /* beta = <r', z> / <r, z> */
      solveur.Dotc(arg, gradr, gradz, alpha, 0);
      solveur.Div(arg, alpha, beta, beta, 1);

      /* p = z + beta * p */
      solveur.BYPX(arg, beta, gradz, gradp, 1);

      normr = solveur.Norm(arg, gradr);
      tmp = normr / normb;

      RAFF_CLOCK_STOP;
      t3 = RAFF_CLOCK_GET;
      solveur.Verbose(arg, t0, t3, tmp, nb_iter);
      t0 = t3;
    }

  solveur.End(arg, tmp, nb_iter, t3, gradx);

  solveur.Free(arg, (void*) gradb);
  solveur.Free(arg, (void*) gradr);
  solveur.Free(arg, (void*) gradp);
  solveur.Free(arg, (void*) gradz);
  solveur.Free(arg, (void*) grad2);
  solveur.Free(arg, (void*) alpha);
  solveur.Free(arg, (void*) beta);
  solveur.Free(arg, (void*) gradx);
  return 0;
}

/*
 ** Section: Function creating threads
 */
void z_grad_thread(SolverMatrix *datacode, SopalinParam *sopaparam)
{
  z_raff_thread(datacode, sopaparam, &z_grad_smp);
}
