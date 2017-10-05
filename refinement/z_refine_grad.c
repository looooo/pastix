/**
 *
 * @file z_refine_grad.c
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "bcsc.h"
#include "z_refine_functions.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_refine
 *
 * z_grad_smp - Refine the solution using conjugate gradian method.
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
void z_grad_smp(pastix_data_t *pastix_data, void *x, void *b)
{
    struct z_solver solveur;
    memset( &solveur, 0, sizeof(struct z_solver) );

    z_Pastix_Solveur(&solveur);

    pastix_bcsc_t     * bcsc      = pastix_data->bcsc;
    pastix_int_t        n         = bcsc->gN;
    Clock               refine_clk;
    pastix_fixdbl_t     t0        = 0;
    pastix_fixdbl_t     t3        = 0;
    pastix_complex64_t  tmp       = 0.0;
    pastix_complex64_t  normr;
    pastix_complex64_t  normb;
    pastix_complex64_t  epsilon   = solveur.Eps(pastix_data);
    pastix_int_t        itermax   = solveur.Itermax(pastix_data);
    pastix_int_t        nb_iter   = 0;

    pastix_complex64_t *gradb = NULL;
    pastix_complex64_t *gradr = NULL;
    pastix_complex64_t *gradp = NULL;
    pastix_complex64_t *gradz = NULL;
    pastix_complex64_t *grad2 = NULL;
    pastix_complex64_t *alpha = NULL;
    pastix_complex64_t *beta  = NULL;
    pastix_complex64_t *gradx = NULL;

    /* Initialize vectors */
    gradb = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    gradr = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    gradp = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    gradz = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    grad2 = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    alpha = (pastix_complex64_t *)solveur.Malloc(    sizeof(pastix_complex64_t));
    beta  = (pastix_complex64_t *)solveur.Malloc(    sizeof(pastix_complex64_t));
    gradx = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));

    clockInit(refine_clk);clockStart(refine_clk);

    solveur.B(b, gradb, n);
    solveur.X(pastix_data, x, gradx);

    /* r=b-ax */
    solveur.bMAx(bcsc, gradb, gradx, gradr);
    normb = solveur.Norm(gradb, n);
    normr = solveur.Norm(gradr, n);
    tmp = normr / normb;

    /* z = M-1 r */
    solveur.Precond(pastix_data, gradr, gradz);

    memcpy(gradp, gradz, n * sizeof( pastix_complex64_t ));

    while (((double)tmp > (double)epsilon) && (nb_iter < itermax))
    {
        clockStop((refine_clk));
        t0 = clockGet();
        nb_iter++;

        /* grad2 = A * p */
        solveur.Ax(bcsc, gradp, grad2);

        /* alpha = <r, z> / <Ap, p> */
        solveur.Dotc(n, gradr, gradz, beta);
        solveur.Dotc(n, grad2, gradp, alpha);
        // solveur.Div(arg, beta, alpha, alpha, 1);
        alpha[0] = beta[0] / alpha[0];

        /* x = x + alpha * p */
        solveur.AXPY(n, 1, alpha, gradx, gradp);

        /* r = r - alpha * A * p */
        solveur.AXPY(n, -1, alpha, gradr, grad2);

        /* z = M-1 * r */
        solveur.Precond(pastix_data, gradr, gradz);

        /* beta = <r', z> / <r, z> */
        solveur.Dotc(n, gradr, gradz, alpha);
        // solveur.Div(arg, alpha, beta, beta, 1);
        beta[0] = alpha[0] / beta[0];

        /* p = z + beta * p */
        solveur.BYPX(n, beta, gradz, gradp);

        normr = solveur.Norm(gradr, n);
        tmp = normr / normb;

        clockStop((refine_clk));
        t3 = clockGet();
        if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot )
            solveur.Verbose(t0, t3, tmp, nb_iter);
        t0 = t3;
    }

    solveur.End(pastix_data, tmp, nb_iter, t3, x, gradx);

    solveur.Free((void*) gradb);
    solveur.Free((void*) gradr);
    solveur.Free((void*) gradp);
    solveur.Free((void*) gradz);
    solveur.Free((void*) grad2);
    solveur.Free((void*) alpha);
    solveur.Free((void*) beta);
    solveur.Free((void*) gradx);
}
