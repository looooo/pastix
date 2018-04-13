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
    pastix_int_t                   im, im1      = 0;
    volatile pastix_int_t          maxits  = 0;
    pastix_complex64_t            *gmHi, *gmH   = NULL;
    pastix_complex64_t            *gmVi, *gmV   = NULL;
    pastix_complex64_t            *gmWi, *gmW   = NULL;
    pastix_complex64_t          *  gmcos        = NULL;
    pastix_complex64_t          *  gmsin        = NULL;
    pastix_complex64_t          *  gmresrs      = NULL;
    volatile pastix_int_t          gmresiters   = 0;
    volatile pastix_int_t          i            = 0;
    pastix_int_t                   j, ii, k, ldw;
    pastix_complex64_t             tmp;
    gmres_t                     *  gmresdata;
    double eps, norm, normb, normx, resid;
    int precond = 1;
    ldw = n;

    if ( bcsc->mtxtype == PastixHermitian ) {
        /* Check if we need dotu for non hermitian matrices (CEA patch) */
        solveur.dot = &z_Pastix_Dotc;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        precond = 0;
        ldw = 0;
    }

    precond = 0;
    ldw = 0;

    /* Get the parameters */
    im     = solveur.Krylov_Space(pastix_data);
    im1    = im + 1;
    maxits = solveur.Itermax(pastix_data);
    eps    = solveur.Eps(pastix_data);

    gmcos     = (pastix_complex64_t *)solveur.Malloc(im  * sizeof(pastix_complex64_t));
    gmsin     = (pastix_complex64_t *)solveur.Malloc(im  * sizeof(pastix_complex64_t));
    gmresrs   = (pastix_complex64_t *)solveur.Malloc(im1 * sizeof(pastix_complex64_t));
    gmresdata = (gmres_t *)solveur.Malloc(1 * sizeof(gmres_t));
    /**
     * H stores the h_{i,j} elements ot the upper hessenberg matrix H (See Alg. 9.5 p 270)
     * V stores the v_{i} vectors
     * W stores the M^{-1} v_{i} vectors to avoid the application of the
     *          preconditioner on the output result (See line 11 of Alg 9.5)
     */
    gmH = (pastix_complex64_t *)solveur.Malloc(im * im1 * sizeof(pastix_complex64_t));
    gmV = (pastix_complex64_t *)solveur.Malloc(n  * im1 * sizeof(pastix_complex64_t));
    if (precond) {
        gmW = (pastix_complex64_t *)solveur.Malloc(n  * im  * sizeof(pastix_complex64_t));
    }
    else {
        gmW = (pastix_complex64_t *)solveur.Malloc(n        * sizeof(pastix_complex64_t));
    }
    memset( gmH, 0, im * im1 * sizeof(pastix_complex64_t*) );

    gmresdata->gmresro = 0.0;
    gmresdata->gmresout_flag = 1;

    normb = solveur.norm( n, b );
    normx = solveur.norm( n, x );

    gmresiters = 0;

    clockInit(refine_clk);
    clockStart(refine_clk);

    /**
     * Algorithm from Iterative Methods for Sparse Linear systems, Y. Saad, Second Ed. p267-273
     *
     * The version implemented is the Right preconditioned algorithm.
     */
    while (gmresdata->gmresout_flag)
    {
        /* Initialize v_{0} and w_{0} */
        gmVi = gmV;
        gmWi = gmW;

        /* Compute r0 = b - A * x */
        solveur.copy( n, b, gmVi );
        if ( normx > 0. ) {
            solveur.spmv( pastix_data, -1., x, 1., gmVi );
        }

        /* Compute beta = ||r0||_f */
        gmresdata->gmresro = solveur.norm( n, gmVi );

        /* If residual is small enough, exit */
        if ( gmresdata->gmresro <= eps )
        {
            gmresdata->gmresout_flag = 0;
            break;
        }

        /* Compute v0 = r0 / beta */
        tmp = (pastix_complex64_t)( 1.0 / gmresdata->gmresro );
        solveur.scal( n, tmp, gmVi );

        gmresrs[0] = (pastix_complex64_t)gmresdata->gmresro;
        gmresdata->gmresin_flag = 1;
        i = -1;
        gmHi = gmH - im1;
        gmWi = gmW - ldw;

        while( gmresdata->gmresin_flag )
        {
            clockStop( refine_clk );
            t0 = clockGet();

            i++;

            /* Set H and W pointers to the beginning of columns i */
            gmHi = gmHi + im1;
            gmWi = gmWi + ldw;

            /* Backup v_{i} into w_{i} for the end */
            solveur.copy( n, gmVi, gmWi );

            /* Compute w_{i} = M^{-1} v_{i} */
            if ( precond ) {
                solveur.spsv( pastix_data, gmWi );
            }

            /* v_{i+1} = A (M^{-1} v_{i}) = A w_{i} */
            gmVi += n;
            solveur.spmv( pastix_data, 1.0, gmWi, 0., gmVi );

            /* Classical Gram-Schmidt */
            for (j=0; j<=i; j++)
            {
                /* Compute h_{j,i} = < v_{i+1}, v_{j} > */
                gmHi[j] = solveur.dot( n, gmVi, gmV + j * n );

                /* Compute v_{i+1} = v_{i+1} - h_{j,i} v_{j} */
                solveur.axpy( n, -1. * gmHi[j],  gmV + j * n, gmVi );
            }

            /* Compute || v_{i+1} ||_f */
            norm = solveur.norm( n, gmVi );
            gmHi[i+1] = norm;

            /* Compute V_{i+1} = v_{i+1} / h_{i+1,i} iff h_{i+1,i} is not too small */
            if ( norm > 1e-50 )
            {
                tmp = (pastix_complex64_t)(1.0 / norm);
                solveur.scal( n, tmp, gmVi );
            }

            /* Apply the previous Givens rotation to the new column (should call LAPACKE_zrot_work())*/
            for (j=0; j<i;j++)
            {
                /*
                 * h_{j,  i} = cos_j * h_{j,  i} +      sin_{j}  * h_{j+1, i}
                 * h_{j+1,i} = cos_j * h_{j+1,i} - conj(sin_{j}) * h_{j,   i}
                 */
                tmp = gmHi[j];
                gmHi[j]   = gmcos[j] * tmp       +      gmsin[j]  * gmHi[j+1];
                gmHi[j+1] = gmcos[j] * gmHi[j+1] - conj(gmsin[j]) * tmp;
            }

            /*
             * Compute the new Givens rotation (zrotg)
             *
             * t = sqrt( h_{i,i}^2 + h_{i+1,i}^2 )
             * cos = h_{i,i}   / t
             * sin = h_{i+1,i} / t
             */
            {
                tmp = csqrt( gmHi[i]   * gmHi[i] +
                             gmHi[i+1] * gmHi[i+1] );

                if ( cabs(tmp) <= eps ) {
                    tmp = (pastix_complex64_t)eps;
                }
                gmcos[i] = gmHi[i]   / tmp;
                gmsin[i] = gmHi[i+1] / tmp;
            }

            /* Update the residuals (See p. 168, eq 6.35) */
            gmresrs[i+1] = -gmsin[i] * gmresrs[i];
            gmresrs[i]   =  gmcos[i] * gmresrs[i];

            /* Apply the last Givens rotation */
            gmHi[i] = gmcos[i] * gmHi[i] + gmsin[i] * gmHi[i+1];

            /* (See p. 169, eq 6.42) */
            gmresdata->gmresro = cabs( gmresrs[i+1] );

            resid = gmresdata->gmresro / normb;
            gmresiters++;
            if ( (i+1 >= im) ||
                 (resid <= eps) ||
                 (gmresiters >= maxits) )
            {
                gmresdata->gmresin_flag = 0;
            }

            clockStop((refine_clk));
            t3 = clockGet();
            if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
                solveur.output_oneiter( t0, t3, resid, gmresiters );
            }
        }

        /* Compute y_m = H_m^{-1} g_m (See p. 169) */
        gmresrs[i] = gmresrs[i] / gmHi[i];
        for (ii=0; ii<i; ii++)
        {
            k = i-ii-1;
            tmp = gmresrs[k];
            for (j=k+1; j<=i; j++)
            {
                tmp = tmp - gmH[j * im1 + k] * gmresrs[j];
            }
            gmresrs[k] = tmp / gmH[k * im1 + k];
        }

        /**
         * TODO: - Replace by Gemv, and add last Precond for right precond
         *       - Add the left preconditionning variant
         * Compute x_m = x_0 + M^{-1} V_m y_m
         *             = x_0 +        W_m y_m
         */
        for (j=0; j<=i;j++)
        {
            if (precond) {
                gmWi = gmW;
            }
            else {
                gmWi = gmV;
            }
            /* X + */
            solveur.axpy( n, gmresrs[j], gmWi + j * n, x );
        }

        if ((resid <= eps) || (gmresiters >= maxits))
        {
            gmresdata->gmresout_flag = 0;
        }
    }

    clockStop( refine_clk );
    t3 = clockGet();

    solveur.End( pastix_data, resid, gmresiters, t3, x, x );

    solveur.Free(gmcos);
    solveur.Free(gmsin);
    solveur.Free(gmresrs);
    solveur.Free(gmresdata);
    solveur.Free(gmH);
    solveur.Free(gmV);
    solveur.Free(gmW);
}
