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
#include "common.h"
#include "bcsc.h"
#include "z_raff_functions.h"
#include "solver.h"

/*
 ** Section: Threads routines
 */

/*
 Function: z_gmres_smp

 Function computing GMRES iterative reffinement.

 Parameters:
 arg - Pointer to a <sopthread_data_t> structure containing
 the <z_Sopalin_Data_t> structure and the thread number ID.
 */

typedef struct gmres_s
{
  volatile pastix_int_t gmresout_flag;     /*+ Flag for GMRES outter loop          +*/
  volatile pastix_int_t gmresin_flag;      /*+ Flag for GMRES inner loop           +*/
  volatile double     gmresro;           /*+ Norm of GMRES residue               +*/
} gmres_t;

void z_gmres_smp(pastix_data_t *pastix_data, void *x, void *b)
{
  struct z_solver solveur = {NULL};
  z_Pastix_Solveur(&solveur);
//   SolverMatrix                 * datacode     = pastix_data->solvmatr;
  SopalinParam                 * sopaparam    = &(pastix_data->sopar);
  pastix_bcsc_t                * bcsc         = pastix_data->bcsc;
  pastix_int_t                   n            = bcsc->gN;
  Clock                          raff_clk;
  double                         t0           = 0.0;
  double                         t3           = 0.0;
  pastix_complex64_t          *  gmrestemp    = NULL;
  volatile pastix_int_t          gmresim      = 0;
  volatile pastix_int_t          gmresmaxits  = 0;
  pastix_complex64_t          *  gmresb       = NULL;
  pastix_complex64_t          ** gmresvv      = NULL;
  pastix_complex64_t          ** gmreshh      = NULL;
  pastix_complex64_t          *  gmresc       = NULL;
  pastix_complex64_t          *  gmress       = NULL;
  pastix_complex64_t          *  gmresrs      = NULL;
  pastix_complex64_t          ** gmresw       = NULL;
  pastix_complex64_t             gmresalpha;
  pastix_complex64_t             gmrest;
  volatile pastix_int_t          gmresiters   = 0;
  pastix_complex64_t          *  gmreswk1;
  pastix_complex64_t          *  gmreswk2     = NULL;
  volatile double                gmreseps     = 0;
  volatile double                gmresnormb;
  volatile pastix_int_t          gmresi1      = 0;
  volatile pastix_int_t          i            = 0;
  pastix_int_t                   j, ii, k;
  pastix_complex64_t             beta;
  pastix_complex64_t          *  gmresx       = NULL;
  gmres_t                     *  gmresdata;
  gmresim     = solveur.Krylov_Space(sopaparam);
  gmresmaxits = solveur.Itermax(sopaparam);
  gmreseps    = solveur.Eps(sopaparam);

  gmrestemp = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
  gmresb    = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
  gmresc    = (pastix_complex64_t *)solveur.Malloc(gmresim * sizeof(pastix_complex64_t));
  gmress    = (pastix_complex64_t *)solveur.Malloc(gmresim * sizeof(pastix_complex64_t));
  gmresrs   = (pastix_complex64_t *)solveur.Malloc((gmresim+1) * sizeof(pastix_complex64_t));
  gmresdata = (gmres_t *)solveur.Malloc(1 * sizeof(gmres_t));
  gmresx    = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
  gmresvv   = (pastix_complex64_t **)solveur.Malloc((gmresim+1) * sizeof(pastix_complex64_t*));
  gmreshh   = (pastix_complex64_t **)solveur.Malloc(gmresim * sizeof(pastix_complex64_t*));
  gmresw    = (pastix_complex64_t **)solveur.Malloc(gmresim * sizeof(pastix_complex64_t*));
  for (i=0; i<gmresim; i++)
    {
      gmresvv[i] = (pastix_complex64_t *)solveur.Malloc(n           * sizeof(pastix_complex64_t));
      gmreshh[i] = (pastix_complex64_t *)solveur.Malloc((gmresim+1) * sizeof(pastix_complex64_t));
      gmresw[i]  = (pastix_complex64_t *)solveur.Malloc(n           * sizeof(pastix_complex64_t));
    }
  gmresvv[gmresim] = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
//   MONO_END(arg);
//   SYNCHRO(arg);

  /* Synchronisations */
//   gmrestemp  = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmrestemp, 0);
//   gmresb     = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmresb,    1);
//   gmresc     = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmresc,    2);
//   gmress     = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmress,    3);
//   gmresrs    = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmresrs,   4);
//   gmresvv    = (pastix_complex64_t **)solveur.Synchro(arg, (void*) gmresvv,   6);
//   gmreshh    = (pastix_complex64_t **)solveur.Synchro(arg, (void*) gmreshh,   7);
//   gmresw     = (pastix_complex64_t **)solveur.Synchro(arg, (void*) gmresw,    8);
//   gmresdata  = (gmres_t*)solveur.Synchro(arg, (void*) gmresdata, 9);
// 
//   gmresnormb = (double)(*((double*)solveur.Synchro(arg, (void*) &gmresnormb, 10)));
//   gmresx     = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmresx,    11);

  gmresdata->gmresro = 0.0;
  gmresdata->gmresout_flag = 1;

  solveur.B(b, gmresb, n);
  gmresnormb = solveur.Norm(gmresb, n);

  solveur.X(pastix_data, x, gmresx);

  gmresalpha = -1.0;
  gmresiters = 0;

//   RAFF_CLOCK_INIT;
  clockInit(raff_clk);clockStart(raff_clk);

  while (gmresdata->gmresout_flag)
    {
      gmreswk2 = gmresvv[0];

      /* gmresvv[0] = b - A * x */
      solveur.bMAx(bcsc, gmresb, gmresx, gmresvv[0]);

      /* ro = vv[0].vv[0] */
      solveur.Dotc_Gmres(n, gmresvv[0], gmresvv[0], &beta, 0);

#if defined(PRECISION_z) || defined(PRECISION_c)
      gmresdata->gmresro = (pastix_complex64_t)csqrt(beta);
#else
      gmresdata->gmresro = (pastix_complex64_t)sqrt(beta);
#endif

      if ((double)cabs((pastix_complex64_t)gmresdata->gmresro) <=
          gmreseps)
        {
          gmresdata->gmresout_flag = 0;
          break;
        }

      gmrest = (pastix_complex64_t)(1.0/gmresdata->gmresro);

      solveur.Scal(n, gmrest, gmresvv[0], 1);

      gmresrs[0] = (pastix_complex64_t)gmresdata->gmresro;
      gmresdata->gmresin_flag = 1;
      i=-1;

      while(gmresdata->gmresin_flag)
        {
//           RAFF_CLOCK_STOP;
          clockStop((raff_clk));
//           t0 = RAFF_CLOCK_GET;
          t0 = clockGet();

          i++;
          gmresi1 = i+1;

          gmreswk1 = gmresvv[i];
          gmreswk2 = gmresw[i];

//           SYNCHRO(arg);
          solveur.Precond(pastix_data, gmreswk1, gmreswk2, 1);

          gmreswk1 = gmresvv[gmresi1];

          /* vv[i1] = A*wk2 */
          solveur.Ax(bcsc, gmreswk2, gmreswk1);

          /* classical gram - schmidt */
          for (j=0; j<=i; j++)
            {
              /* vv[j]*vv[i1] */
              solveur.Dotc_Gmres(n, gmresvv[gmresi1], gmresvv[j], &beta, 0);

              gmreshh[i][j] = (pastix_complex64_t)beta;
            }

//           SYNCHRO(arg);

          for (j=0;j<=i;j++)
            {
              gmresalpha = -gmreshh[i][j];
              solveur.AXPY(n, 1.0, &gmresalpha, gmresvv[gmresi1], gmresvv[j], 0);
            }

//           SYNCHRO(arg);
          solveur.Dotc_Gmres(n, gmresvv[gmresi1], gmresvv[gmresi1], &beta, 0);

#if defined(PRECISION_z) || defined(PRECISION_c)
      gmrest = (pastix_complex64_t)csqrt(beta);
#else
      gmrest = (pastix_complex64_t)sqrt(beta);
#endif

          gmreshh[i][gmresi1] = gmrest;

          if (cabs(gmrest) > 10e-50)
            {
//               gmrest = fun / gmrest;
              gmrest = 1.0 / gmrest;
              solveur.Scal(n, gmrest, gmresvv[gmresi1], 0);
            }

//           SYNCHRO(arg);
//           MONO_BEGIN(arg);

          if (i != 0)
            {
              for (j=1; j<=i;j++)
                {
                  gmrest = gmreshh[i][j-1];
#ifdef TYPE_COMPLEX
                  gmreshh[i][j-1] = (pastix_complex64_t)conj(gmresc[j-1])*gmrest +
                    (pastix_complex64_t)conj(gmress[j-1])*gmreshh[i][j];
#else /* TYPE_COMPLEX */
                  gmreshh[i][j-1] =  gmresc[j-1]*gmrest +
                    gmress[j-1]*gmreshh[i][j];
#endif /* TYPE_COMPLEX */
                  gmreshh[i][j]   = -gmress[j-1]*gmrest +
                    gmresc[j-1]*gmreshh[i][j];
                }
            }
#ifdef TYPE_COMPLEX
          gmrest = (pastix_complex64_t)csqrt(cabs(gmreshh[i][i]*gmreshh[i][i])+
                                     gmreshh[i][gmresi1]*gmreshh[i][gmresi1]);
#else
          gmrest = (pastix_complex64_t)sqrt(gmreshh[i][i]*gmreshh[i][i]+
                                    gmreshh[i][gmresi1]*gmreshh[i][gmresi1]);
#endif
          if (cabs(gmrest) <= gmreseps)
            gmrest = (pastix_complex64_t)gmreseps;

          gmresc[i] = gmreshh[i][i]/gmrest;
          gmress[i] = gmreshh[i][gmresi1]/gmrest;
          gmresrs[gmresi1] = -gmress[i]*gmresrs[i];

#ifdef TYPE_COMPLEX
          gmresrs[i] = (pastix_complex64_t)conj(gmresc[i])*gmresrs[i];
          gmreshh[i][i] = (pastix_complex64_t)conj(gmresc[i])*gmreshh[i][i] +
          gmress[i]*gmreshh[i][gmresi1];
#else
          gmresrs[i] = gmresc[i]*gmresrs[i];
          gmreshh[i][i] = gmresc[i]*gmreshh[i][i] +
          gmress[i]*gmreshh[i][gmresi1];
#endif
          gmresdata->gmresro = cabs(gmresrs[gmresi1]);

//           MONO_END(arg);

          gmresiters++;

//           MONO_BEGIN(arg);
          if ((i+1 >= gmresim) || (gmresdata->gmresro/gmresnormb <= gmreseps) || (gmresiters >= gmresmaxits))
            {
              gmresdata->gmresin_flag = 0;
            }
//           MONO_END(arg);

//           RAFF_CLOCK_STOP;
          clockStop((raff_clk));
//           t3 = RAFF_CLOCK_GET;
          t3 = clockGet();
          solveur.Verbose(t0, t3, gmresdata->gmresro/gmresnormb, gmresiters);
//           SYNCHRO(arg);
        }

//       MONO_BEGIN(arg);

      gmresrs[i] = gmresrs[i]/gmreshh[i][i];
      for (ii=2; ii<=i+1; ii++)
        {
          k = i-ii+1;
          gmrest = gmresrs[k];
          for (j=k+1; j<=i; j++)
            {
              gmrest = gmrest - gmreshh[j][k]*gmresrs[j];
            }
          gmresrs[k] = gmrest/gmreshh[k][k];
        }

//       MONO_END(arg);
//       SYNCHRO(arg);

      for (j=0; j<=i;j++)
        {
          gmrest = gmresrs[j];
          solveur.AXPY(n, 1.0, &gmrest, gmresx, gmresw[j], 0);
        }
//       SYNCHRO(arg);

      if ((gmresdata->gmresro/gmresnormb<= gmreseps) || (gmresiters >= gmresmaxits))
        {
          gmresdata->gmresout_flag = 0;
        }
    }

//   RAFF_CLOCK_STOP;
  clockStop((raff_clk));
//   t3 = RAFF_CLOCK_GET;
  t3 = clockGet();

  solveur.End(sopaparam, gmresdata->gmresro/gmresnormb, gmresiters, t3, x, gmresx);

  solveur.Free((void*) gmrestemp);
  solveur.Free((void*) gmresb);
  solveur.Free((void*) gmresc);
  solveur.Free((void*) gmress);
  solveur.Free((void*) gmresrs);
  solveur.Free((void*) gmresdata);
  solveur.Free((void*) gmresx);

//   MONO_BEGIN(arg);
  for (i=0; i<gmresim; i++)
    {
      solveur.Free(gmresvv[i]);
      solveur.Free(gmreshh[i]);
      solveur.Free(gmresw[i]);
    }

  solveur.Free(gmresvv[gmresim]);

  solveur.Free(gmresvv);
  solveur.Free(gmreshh);
  solveur.Free(gmresw);
//   MONO_END(arg);

//   return 0;
}

/*
** Section: Function creating threads
*/
// void gmres_thread(pastix_data_t *pastix_data, void *x, void *b)
// {
// //   z_raff_thread(bcsc, sopaparam, &z_gmres_smp);
//     int i;
//     double *bptr = (double*)b;
//     for(i=0; i<pastix_data->bcsc->gN*2; i+=2, bptr+=2)
//         fprintf(stdout, "++     b[0] = %g +I%g\n", bptr[0],bptr[1]);
//     z_gmres_smp(pastix_data, x, b);
// }
