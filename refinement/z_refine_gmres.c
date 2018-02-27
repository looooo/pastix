/**
 *
 * @file z_refine_gmres.c
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Theophile Terraz
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "bcsc.h"
#include "z_refine_functions.h"

/**
 *
 * @brief GMRES parameters
 *
 * This structure describes gmres parameters.
 *
 */
typedef struct gmres_s
{
    volatile pastix_int_t gmresout_flag;     /*+ Flag for GMRES outter loop          +*/
    volatile pastix_int_t gmresin_flag;      /*+ Flag for GMRES inner loop           +*/
    volatile double       gmresro;           /*+ Norm of GMRES residue               +*/
} gmres_t;

/**
 *******************************************************************************
 *
 * @ingroup pastix_refine
 *
 * z_gmres_smp - Function computing GMRES iterative refinement.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[out] x
 *          The solution vector.
 *
 * @param[in] b
 *          The right hand side member (only one).
 *
 *******************************************************************************/
void z_gmres_smp(pastix_data_t *pastix_data, void *x, void *b)
{
    struct z_solver solveur;
    memset( &solveur, 0, sizeof(struct z_solver) );

    z_Pastix_Solveur(&solveur);

    pastix_bcsc_t                * bcsc         = pastix_data->bcsc;
    pastix_int_t                   n            = bcsc->gN;
    Clock                          refine_clk;
    pastix_fixdbl_t                t0           = 0.0;
    pastix_fixdbl_t                t3           = 0.0;
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

    if(bcsc->mtxtype == PastixHermitian)
        solveur.Dotc = &z_Pastix_Dotu;
    gmresim     = solveur.Krylov_Space(pastix_data);
    gmresmaxits = solveur.Itermax(pastix_data);
    gmreseps    = solveur.Eps(pastix_data);

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
    gmresdata->gmresro = 0.0;
    gmresdata->gmresout_flag = 1;

    solveur.B(b, gmresb, n);
    gmresnormb = solveur.Norm(gmresb, n);

    solveur.X(pastix_data, x, gmresx);

    gmresalpha = -1.0;
    gmresiters = 0;

    clockInit(refine_clk);clockStart(refine_clk);

    while (gmresdata->gmresout_flag)
    {
        gmreswk2 = gmresvv[0];

        /* gmresvv[0] = b - A * x */
        solveur.bMAx(bcsc, gmresb, gmresx, gmresvv[0]);

        /* ro = vv[0].vv[0] */
        solveur.Dotc(n, gmresvv[0], gmresvv[0], &beta);

        gmresdata->gmresro = (pastix_complex64_t)csqrt(beta);

        if ((double)cabs((pastix_complex64_t)gmresdata->gmresro) <=
            gmreseps)
        {
            gmresdata->gmresout_flag = 0;
            break;
        }

        gmrest = (pastix_complex64_t)(1.0/gmresdata->gmresro);

        solveur.Scal(n, gmrest, gmresvv[0]);

        gmresrs[0] = (pastix_complex64_t)gmresdata->gmresro;
        gmresdata->gmresin_flag = 1;
        i=-1;

        while(gmresdata->gmresin_flag)
        {
            clockStop((refine_clk));
            t0 = clockGet();

            i++;
            gmresi1 = i+1;

            gmreswk1 = gmresvv[i];
            gmreswk2 = gmresw[i];

            solveur.Precond(pastix_data, gmreswk1, gmreswk2);

            gmreswk1 = gmresvv[gmresi1];

            /* vv[i1] = A*wk2 */
            solveur.Ax(bcsc, gmreswk2, gmreswk1);

            /* classical gram - schmidt */
            for (j=0; j<=i; j++)
            {
                /* vv[j]*vv[i1] */
                solveur.Dotc(n, gmresvv[gmresi1], gmresvv[j], &beta);

                gmreshh[i][j] = (pastix_complex64_t)beta;
            }

            for (j=0;j<=i;j++)
            {
                gmresalpha = -gmreshh[i][j];
                solveur.AXPY(n, 1.0, &gmresalpha, gmresvv[gmresi1], gmresvv[j]);
            }

            solveur.Dotc(n, gmresvv[gmresi1], gmresvv[gmresi1], &beta);

            gmrest = (pastix_complex64_t)csqrt(beta);

            gmreshh[i][gmresi1] = gmrest;

            if (cabs(gmrest) > 10e-50)
            {
                //               gmrest = fun / gmrest;
                gmrest = 1.0 / gmrest;
                solveur.Scal(n, gmrest, gmresvv[gmresi1]);
            }

            if (i != 0)
            {
                for (j=1; j<=i;j++)
                {
                    gmrest = gmreshh[i][j-1];
                    gmreshh[i][j-1] = (pastix_complex64_t)conj(gmresc[j-1])*gmrest +
                        (pastix_complex64_t)conj(gmress[j-1])*gmreshh[i][j];
                    gmreshh[i][j]   = -gmress[j-1]*gmrest +
                        gmresc[j-1]*gmreshh[i][j];
                }
            }
            gmrest = (pastix_complex64_t)csqrt(cabs(gmreshh[i][i]*gmreshh[i][i])+
                                               gmreshh[i][gmresi1]*gmreshh[i][gmresi1]);
            if (cabs(gmrest) <= gmreseps)
                gmrest = (pastix_complex64_t)gmreseps;

            gmresc[i] = gmreshh[i][i]/gmrest;
            gmress[i] = gmreshh[i][gmresi1]/gmrest;
            gmresrs[gmresi1] = -gmress[i]*gmresrs[i];

            gmresrs[i] = (pastix_complex64_t)conj(gmresc[i])*gmresrs[i];
            gmreshh[i][i] = (pastix_complex64_t)conj(gmresc[i])*gmreshh[i][i] +
                gmress[i]*gmreshh[i][gmresi1];
            gmresdata->gmresro = cabs(gmresrs[gmresi1]);

            gmresiters++;

            if ((i+1 >= gmresim) || (gmresdata->gmresro/gmresnormb <= gmreseps) || (gmresiters >= gmresmaxits))
            {
                gmresdata->gmresin_flag = 0;
            }

            clockStop((refine_clk));
            t3 = clockGet();
            if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot )
                solveur.Verbose(t0, t3, gmresdata->gmresro/gmresnormb, gmresiters);
        }

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

        for (j=0; j<=i;j++)
        {
            gmrest = gmresrs[j];
            solveur.AXPY(n, 1.0, &gmrest, gmresx, gmresw[j]);
        }

        if ((gmresdata->gmresro/gmresnormb<= gmreseps) || (gmresiters >= gmresmaxits))
        {
            gmresdata->gmresout_flag = 0;
        }
    }

    clockStop((refine_clk));
    t3 = clockGet();

    solveur.End(pastix_data, gmresdata->gmresro/gmresnormb, gmresiters, t3, x, gmresx);

    solveur.Free((void*) gmresb);
    solveur.Free((void*) gmresc);
    solveur.Free((void*) gmress);
    solveur.Free((void*) gmresrs);
    solveur.Free((void*) gmresdata);
    solveur.Free((void*) gmresx);

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
}
