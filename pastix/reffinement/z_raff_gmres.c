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
#include "bcsc.h"
#include "z_raff_functions.h"

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

void* z_gmres_smp(void *arg)
{
  struct z_solver solveur = {NULL};
  z_Pastix_Solveur(&solveur);

  pastix_int_t               n            = solveur.N(arg);
  Clock             raff_clk;
  double            t0           = 0;
  double            t3           = 0;
  pastix_complex64_t          *  gmrestemp    = NULL;
  volatile pastix_int_t      gmresim      = 0;
  volatile pastix_int_t      gmresmaxits  = 0;
  pastix_complex64_t          *  gmresb       = NULL;
  pastix_complex64_t          ** gmresvv      = NULL;
  pastix_complex64_t          ** gmreshh      = NULL;
  pastix_complex64_t          *  gmresc       = NULL;
  pastix_complex64_t          *  gmress       = NULL;
  pastix_complex64_t          *  gmresrs      = NULL;
  pastix_complex64_t          ** gmresw       = NULL;
  pastix_complex64_t             gmresalpha;
  pastix_complex64_t             gmrest;
  volatile pastix_int_t      gmresiters   = 0;
  pastix_complex64_t          *  gmreswk1;
  pastix_complex64_t          *  gmreswk2     = NULL;
  volatile double   gmreseps     = 0;
  volatile double   gmresnormb;
  volatile pastix_int_t      gmresi1      = 0;
  volatile pastix_int_t      i = 0;
  pastix_int_t               j, ii, k;
  pastix_complex64_t             beta;
  pastix_complex64_t          *  gmresx       = NULL;
  gmres_t        *  gmresdata;
  gmresim     = solveur.Krylov_Space(arg);
  gmresmaxits = solveur.Itermax(arg);
  gmreseps    = solveur.Eps(arg);

  gmrestemp = solveur.Malloc(arg, n           * sizeof(pastix_complex64_t));
  gmresb    = solveur.Malloc(arg, n           * sizeof(pastix_complex64_t));
  gmresc    = solveur.Malloc(arg, gmresim     * sizeof(pastix_complex64_t));
  gmress    = solveur.Malloc(arg, gmresim     * sizeof(pastix_complex64_t));
  gmresrs   = solveur.Malloc(arg, (gmresim+1) * sizeof(pastix_complex64_t));
  gmresdata = solveur.Malloc(arg, 1           * sizeof(gmres_t));
  gmresx    = solveur.Malloc(arg, n           * sizeof(pastix_complex64_t));

//   MONO_BEGIN(arg);
  gmresvv = solveur.Malloc(arg, (gmresim+1) * sizeof(pastix_complex64_t*));
  gmreshh = solveur.Malloc(arg, gmresim     * sizeof(pastix_complex64_t*));
  gmresw  = solveur.Malloc(arg, gmresim     * sizeof(pastix_complex64_t*));
  for (i=0; i<gmresim; i++)
    {
      gmresvv[i] = solveur.Malloc(arg, n           * sizeof(pastix_complex64_t));
      gmreshh[i] = solveur.Malloc(arg, (gmresim+1) * sizeof(pastix_complex64_t));
      gmresw[i]  = solveur.Malloc(arg, n           * sizeof(pastix_complex64_t));
    }
  gmresvv[gmresim] = solveur.Malloc(arg, n * sizeof(pastix_complex64_t));
//   MONO_END(arg);
//   SYNCHRO(arg);

  /* Synchronisations */
  gmrestemp  = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmrestemp, 0);
  gmresb     = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmresb,    1);
  gmresc     = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmresc,    2);
  gmress     = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmress,    3);
  gmresrs    = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmresrs,   4);
  gmresvv    = (pastix_complex64_t **)solveur.Synchro(arg, (void*) gmresvv,   6);
  gmreshh    = (pastix_complex64_t **)solveur.Synchro(arg, (void*) gmreshh,   7);
  gmresw     = (pastix_complex64_t **)solveur.Synchro(arg, (void*) gmresw,    8);
  gmresdata  = (gmres_t*)solveur.Synchro(arg, (void*) gmresdata, 9);

  gmresnormb = (double)(*((double*)solveur.Synchro(arg, (void*) &gmresnormb, 10)));
  gmresx     = (pastix_complex64_t * )solveur.Synchro(arg, (void*) gmresx,    11);

  gmresdata->gmresro = 0.0;
  gmresdata->gmresout_flag = 1;

  solveur.B(arg, gmresb);
  gmresnormb = solveur.Norm(arg, gmresb);

  solveur.X(arg, gmresx);

  gmresalpha = -1.0;
  gmresiters = 0;

//   RAFF_CLOCK_INIT;

  while (gmresdata->gmresout_flag)
    {
      gmreswk2 = gmresvv[0];

      /* gmresvv[0] = b - A * x */
      solveur.bMAx(arg, gmresb, gmresx, gmresvv[0]);

      /* ro = vv[0].vv[0] */
      solveur.Dotc_Gmres(arg, gmresvv[0], gmresvv[0], &beta, 0);

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

      solveur.Scal(arg, gmrest, gmresvv[0], 1);

      gmresrs[0] = (pastix_complex64_t)gmresdata->gmresro;
      gmresdata->gmresin_flag = 1;
      i=-1;

      while(gmresdata->gmresin_flag)
        {
//           RAFF_CLOCK_STOP;
//           t0 = RAFF_CLOCK_GET;

          i++;
          gmresi1 = i+1;

          gmreswk1 = gmresvv[i];
          gmreswk2 = gmresw[i];

          SYNCHRO(arg);
          solveur.Precond(arg, gmreswk1, gmreswk2, 1);

          gmreswk1 = gmresvv[gmresi1];

          /* vv[i1] = A*wk2 */
          solveur.Ax(arg, gmreswk2, gmreswk1);

          /* classical gram - schmidt */
          for (j=0; j<=i; j++)
            {
              /* vv[j]*vv[i1] */
              solveur.Dotc_Gmres(arg,gmresvv[gmresi1], gmresvv[j], &beta, 0);

              gmreshh[i][j] = (pastix_complex64_t)beta;
            }

//           SYNCHRO(arg);

          for (j=0;j<=i;j++)
            {
              gmresalpha = -gmreshh[i][j];
              solveur.AXPY(arg, 1.0, &gmresalpha, gmresvv[gmresi1], gmresvv[j], 0);
            }

//           SYNCHRO(arg);
          solveur.Dotc_Gmres(arg, gmresvv[gmresi1], gmresvv[gmresi1], &beta, 0);

#if defined(PRECISION_z) || defined(PRECISION_c)
      gmrest = (pastix_complex64_t)csqrt(beta);
#else
      gmrest = (pastix_complex64_t)sqrt(beta);
#endif

          gmreshh[i][gmresi1] = gmrest;

          if (cabs(gmrest) > 10e-50)
            {
              gmrest = fun / gmrest;
              solveur.Scal(arg, gmrest, gmresvv[gmresi1], 0);
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

          MONO_BEGIN(arg);
          if ((i+1 >= gmresim) || (gmresdata->gmresro/gmresnormb <= gmreseps) || (gmresiters >= gmresmaxits))
            {
              gmresdata->gmresin_flag = 0;
            }
//           MONO_END(arg);

          RAFF_CLOCK_STOP;
          t3 = RAFF_CLOCK_GET;
          solveur.Verbose(arg, t0, t3, gmresdata->gmresro/gmresnormb, gmresiters);
          SYNCHRO(arg);
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
          solveur.AXPY(arg, 1.0, &gmrest, gmresx, gmresw[j], 0);
        }
//       SYNCHRO(arg);

      if ((gmresdata->gmresro/gmresnormb<= gmreseps) || (gmresiters >= gmresmaxits))
        {
          gmresdata->gmresout_flag = 0;
        }
    }

//   RAFF_CLOCK_STOP;
//   t3 = RAFF_CLOCK_GET;

  solveur.End(arg, gmresdata->gmresro/gmresnormb, gmresiters, t3, gmresx);

  solveur.Free(arg, (void*) gmrestemp);
  solveur.Free(arg, (void*) gmresb);
  solveur.Free(arg, (void*) gmresc);
  solveur.Free(arg, (void*) gmress);
  solveur.Free(arg, (void*) gmresrs);
  solveur.Free(arg, (void*) gmresdata);
  solveur.Free(arg, (void*) gmresx);

//   MONO_BEGIN(arg);
  for (i=0; i<gmresim; i++)
    {
      solveur.Free(arg, gmresvv[i]);
      solveur.Free(arg, gmreshh[i]);
      solveur.Free(arg, gmresw[i]);
    }

  solveur.Free(arg, gmresvv[gmresim]);

  solveur.Free(arg, gmresvv);
  solveur.Free(arg, gmreshh);
  solveur.Free(arg, gmresw);
//   MONO_END(arg);

  return 0;
}

/*
** Section: Function creating threads
*/
void gmres_thread(pastix_bcsc_t *bcsc, SopalinParam *sopaparam)
{
  z_raff_thread(bcsc, sopaparam, &z_gmres_smp);
}
