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
    pastix_int_t                   n            = bcsc->n;
    Clock                          refine_clk;
    pastix_fixdbl_t                t0           = 0.0;
    pastix_fixdbl_t                t3           = 0.0;
    volatile pastix_int_t          gmresim      = 0;
    volatile pastix_int_t          gmresmaxits  = 0;
    pastix_complex64_t          ** gmresvv      = NULL;
    pastix_complex64_t          ** gmreshh      = NULL;
    pastix_complex64_t          *  gmcos        = NULL;
    pastix_complex64_t          *  gmsin        = NULL;
    pastix_complex64_t          *  gmresrs      = NULL;
    pastix_complex64_t          ** gmresw       = NULL;
    pastix_complex64_t             gmresalpha;
    pastix_complex64_t             gmrest;
    volatile pastix_int_t          gmresiters   = 0;
    pastix_complex64_t          *  gmreswk1;
    pastix_complex64_t          *  gmreswk2     = NULL;
    volatile double                gmreseps     = 0;
    volatile double                gmresnormb;
    volatile pastix_int_t          i            = 0;
    pastix_int_t                   j, ii, k;
    pastix_complex64_t             beta;
    pastix_complex64_t          *  gmresx       = NULL;
    gmres_t                     *  gmresdata;
    double norm;

    if ( bcsc->mtxtype == PastixHermitian ) {
        /* Check if we need dotu for non hermitian matrices (CEA patch) */
        solveur.Dotc = &z_Pastix_Dotc;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        solveur.Precond = NULL;
    }
    solveur.Precond = NULL;

    /* Get the parameters */
    gmresim     = solveur.Krylov_Space(pastix_data);
    gmresmaxits = solveur.Itermax(pastix_data);
    gmreseps    = solveur.Eps(pastix_data);

    gmcos     = (pastix_complex64_t *)solveur.Malloc(gmresim * sizeof(pastix_complex64_t));
    gmsin     = (pastix_complex64_t *)solveur.Malloc(gmresim * sizeof(pastix_complex64_t));
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

    gmresnormb = solveur.Norm( n, b );

    solveur.X( pastix_data, x, gmresx );

    gmresalpha = -1.0;
    gmresiters = 0;

    clockInit(refine_clk);clockStart(refine_clk);

    /**
     * Algorithm from Iterative Methods for Sparse Linear systems, Y. Saad, Second Ed. p267-273
     *
     * The version implemented is the Right preconditionned algorithm.
     */
    while (gmresdata->gmresout_flag)
    {
        /* Compute r0 = b - A * x */
        solveur.bMAx( bcsc, b, gmresx, gmresvv[0] );

        /* Compute beta = ||r0||_f */
        gmresdata->gmresro = solveur.Norm( n, gmresvv[0] );

        /* If residual is small enough, exit */
        if ( gmresdata->gmresro <= gmreseps )
        {
            gmresdata->gmresout_flag = 0;
            break;
        }

        /* Compute v0 = r0 / beta */
        gmrest = (pastix_complex64_t)( 1.0 / gmresdata->gmresro );
        solveur.Scal( n, gmrest, gmresvv[0] );

        gmresrs[0] = (pastix_complex64_t)gmresdata->gmresro;
        gmresdata->gmresin_flag = 1;
        i = -1;

        while( gmresdata->gmresin_flag )
        {
            clockStop( refine_clk );
            t0 = clockGet();

            i++;
            gmreswk1 = gmresvv[i+1];
            gmreswk2 = gmresw[i];

            /* Compute w2 = M^{-1} v_{i} */
            solveur.B( gmresvv[i], gmreswk2, n );
            if ( solveur.Precond != NULL ) {
                solveur.Precond( pastix_data, gmreswk2 );
            }

            /* w = A (M^{-1} v_{i}) = A w2 */
            solveur.Ax( bcsc, gmreswk2, gmreswk1 );

            /* Classical Gram-Schmidt */
            for (j=0; j<=i; j++)
            {
                /* Compute h_{i,j} = < w, v_{i} > */
                solveur.Dotc( n, gmreswk1, gmresvv[j], &beta );

                gmreshh[i][j] = (pastix_complex64_t)beta;
                gmresalpha = -1. * gmreshh[i][j];

                /* Compute w = w - h_{i,j} v_{i} */
                solveur.AXPY( n, gmresalpha, gmresvv[j], gmreswk1 );
            }

            /* Compute || w ||_f */
            norm = solveur.Norm( n, gmreswk1 );
            gmreshh[i][i+1] = norm;

            /* Compute V_{i+1} = w / h_{i,i+1} iff h_{i,i+1} is not too small */
            if ( norm > 1e-50 )
            {
                gmrest = (pastix_complex64_t)(1.0 / norm);
                solveur.Scal( n, gmrest, gmreswk1 );
            }

            /* Apply the previous Givens rotation to the new column (should call LAPACKE_zrot_work())*/
            if (i != 0)
            {
                for (j=0; j<i;j++)
                {
                    /*
                     * h_{j,  i} = cos_j * h_{j,  i} +      sin_{j}  * h_{j+1, i}
                     * h_{j+1,i} = cos_j * h_{j+1,i} - conj(sin_{j}) * h_{j,   i}
                     */
                    gmrest = gmreshh[i][j];
                    gmreshh[i][j]   = gmcos[j] * gmrest          +      gmsin[j]  * gmreshh[i][j+1];
                    gmreshh[i][j+1] = gmcos[j] * gmreshh[i][j+1] - conj(gmsin[j]) * gmrest;
                }
            }

            /*
             * Compute the new Givens rotation (zrotg)
             *
             * t = sqrt( h_{i,i}^2 + h_{i+1,i}^2 )
             * cos = h_{i,i}   / t
             * sin = h_{i+1,i} / t
             */
            gmrest = csqrt( gmreshh[i][i]   * gmreshh[i][i] +
                            gmreshh[i][i+1] * gmreshh[i][i+1] );

            if ( cabs(gmrest) <= gmreseps ) {
                gmrest = (pastix_complex64_t)gmreseps;
            }
            gmcos[i] = gmreshh[i][i]   / gmrest;
            gmsin[i] = gmreshh[i][i+1] / gmrest;

            /* Update the residuals (See p. 168, eq 6.35) */
            gmresrs[i+1] = -gmsin[i] * gmresrs[i];
            gmresrs[i]   =  gmcos[i] * gmresrs[i];

            /* Apply the last Givens rotation */
            gmreshh[i][i] = gmcos[i] * gmreshh[i][i] + gmsin[i] * gmreshh[i][i+1];

            /* (See p. 169, eq 6.42) */
            gmresdata->gmresro = cabs( gmresrs[i+1] );

            gmresiters++;
            if ( (i+1 >= gmresim) ||
                 ((gmresdata->gmresro / gmresnormb) <= gmreseps) ||
                 (gmresiters >= gmresmaxits) )
            {
                gmresdata->gmresin_flag = 0;
            }

            clockStop((refine_clk));
            t3 = clockGet();
            if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
                solveur.Verbose( t0, t3, gmresdata->gmresro / gmresnormb, gmresiters );
            }
        }

        gmresrs[i] = gmresrs[i] / gmreshh[i][i];
        for (ii=0; ii<i; ii++)
        {
            k = i-ii-1;
            gmrest = gmresrs[k];
            for (j=k+1; j<=i; j++)
            {
                gmrest = gmrest - gmreshh[j][k] * gmresrs[j];
            }
            gmresrs[k] = gmrest / gmreshh[k][k];
        }

        /**
         * TODO: - Replace by Gemv, and add last Precond for right precond
         *       - Add the left preconditionning variant
         */
        for (j=0; j<=i;j++)
        {
            gmrest = gmresrs[j];
            solveur.AXPY( n, gmrest, gmresw[j], gmresx );
        }

        if ((gmresdata->gmresro/gmresnormb<= gmreseps) || (gmresiters >= gmresmaxits))
        {
            gmresdata->gmresout_flag = 0;
        }
    }

    clockStop((refine_clk));
    t3 = clockGet();

    solveur.End(pastix_data, gmresdata->gmresro / gmresnormb, gmresiters, t3, x, gmresx);

    solveur.Free((void*) gmcos);
    solveur.Free((void*) gmsin);
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
