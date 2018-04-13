/**
 *
 * @file z_refine_bicgstab.c
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
 * z_bicgstab_smp - Function computing bicgstab iterative refinement.
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
void z_bicgstab_smp (pastix_data_t *pastix_data, void *x, void *b)
{
    struct z_solver     solver;
    pastix_int_t        n;
    Clock               refine_clk;
    pastix_fixdbl_t     t0      = 0;
    pastix_fixdbl_t     t3      = 0;
    int                 itermax;
    int                 nb_iter = 0;
    int                 precond = 1;
    pastix_complex64_t *gradr; /* Current solution */
    pastix_complex64_t *gradr2;
    pastix_complex64_t *gradp;
    pastix_complex64_t *grady;
    pastix_complex64_t *gradv;
    pastix_complex64_t *grads;
    pastix_complex64_t *gradz;
    pastix_complex64_t *gradt;
    pastix_complex64_t *grad2;
    pastix_complex64_t *grad3;
    pastix_complex64_t  v1, v2, w;
    double normb, normx, normr, alpha, beta;
    double resid_b, eps;

    memset( &solver, 0, sizeof(struct z_solver) );
    z_Pastix_Solver(&solver);

    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        precond = 0;
    }

    n       = solver.getN(    pastix_data );
    itermax = solver.getImax( pastix_data );
    eps     = solver.getEps(  pastix_data );

    gradr  = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    gradr2 = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    gradp  = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    grady  = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    gradv  = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    grads  = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    gradz  = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    gradt  = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    grad2  = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    grad3  = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));

    clockInit(refine_clk);clockStart(refine_clk);

    normb = solver.norm( n, b );
    normx = solver.norm( n, x );

    /* r = b - Ax */
    solver.copy( n, b, gradr );
    if ( normx > 0. ) {
        solver.spmv( pastix_data, -1., x, 1., gradr );
    }
    normr = solver.norm( n, gradr );

    /* r2 = r */
    solver.copy( n, gradr, gradr2 );
    /* p = r */
    solver.copy( n, gradr, gradp );

    /* resid_b = ||r|| / ||b|| */
    resid_b = normr / normb;

    while ((resid_b > eps) && (nb_iter < itermax))
    {
        clockStop((refine_clk));
        t0 = clockGet();
        nb_iter++;

        /* y = M-1 * p */
        solver.copy( n, gradp, grady );
        if ( precond ) {
            solver.spsv( pastix_data, grady );
        }

        /* v = Ay */
        solver.spmv( pastix_data, 1.0, grady, 0., gradv );

        /* alpha = (r, r2) / (v, r2) */
        alpha = solver.dot( n, gradv, gradr2 );
        beta  = solver.dot( n, gradr, gradr2 );
        alpha = beta / alpha;

        /* s = r - alpha * v */
        solver.copy( n, gradr, grads );
        solver.axpy( n, -alpha, gradv, grads );

        /* z = M^{-1} s */
        solver.copy( n, grads, gradz );
        if ( precond ) {
            solver.spsv( pastix_data, gradz );
        }

        /* t = Az */
        solver.spmv( pastix_data, 1.0, gradz, 0., gradt );

        /* w = (M-1t, M-1s) / (M-1t, M-1t) */
        /* grad2 = M-1t */
        solver.copy( n, gradt, grad2 );
        if ( precond ) {
            solver.spsv( pastix_data, grad2 );
        }

        /* v1 = (M-1t, M-1s) */
        /* v2 = (M-1t, M-1t) */
        v1 = solver.dot( n, grad2, gradz );
        v2 = solver.dot( n, grad2, grad2 );
        w = v1 / v2;

        /* x = x + alpha * y + w * z */
        /* x = x + alpha * y */
        solver.axpy( n, alpha, grady, x );

        /* x = x + w * z */
        solver.axpy( n, w, gradz, x );

        /* r = s - w * t*/
        solver.copy( n, grads, gradr );
        solver.axpy( n, -w, gradt, gradr );

        /* beta = (r', r2) / (r, r2) * (alpha / w) */
        /* v1 = (r', r2) */
        v1 = solver.dot( n, gradr, gradr2 );
        v2 = alpha / w;

        beta = v1 / beta;
        beta = beta * v2;

        /* p = r + beta * (p - w * v) */
        /* p = p - w * v */
        solver.axpy( n, -w, gradv, gradp );

        /* p = r + beta * p */
        solver.scal( n, beta, gradp );
        solver.axpy( n, 1., gradr, gradp );

        normr = solver.norm( n, gradr );
        resid_b = normr / normb;

        clockStop((refine_clk));
        t3 = clockGet();
        if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
            solver.output_oneiter( t0, t3, resid_b, nb_iter );
        }
    }

    solver.output_final(pastix_data, resid_b, nb_iter, t3, x, x);

    solver.free((void*) gradr);
    solver.free((void*) gradr2);
    solver.free((void*) gradp);
    solver.free((void*) grady);
    solver.free((void*) gradv);
    solver.free((void*) grads);
    solver.free((void*) gradz);
    solver.free((void*) gradt);
    solver.free((void*) grad2);
    solver.free((void*) grad3);
}
