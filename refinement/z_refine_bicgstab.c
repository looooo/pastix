/**
 *
 * @file z_refine_bicgstab.c
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @author Gregoire Pichon
 * @author Vincent Bridonneau
 * @date 2023-07-21
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "bcsc/bcsc.h"
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
 * @param[out] xp
 *          The solution vector.
 *
 * @param[in] bp
 *          The right hand side member (only one).
 *
 *******************************************************************************
 *
 * @return Number of iterations
 *
 *******************************************************************************/
pastix_int_t
z_bicgstab_smp( pastix_data_t *pastix_data,
                pastix_rhs_t   xp,
                pastix_rhs_t   bp )
{
    struct z_solver     solver;
    pastix_int_t        n;
    Clock               refine_clk;
    pastix_fixdbl_t     t0      = 0;
    pastix_fixdbl_t     t3      = 0;
    int                 itermax;
    int                 nb_iter = 0;
    int                 precond = 1;
    pastix_complex64_t *x = (pastix_complex64_t*)(xp->b);
    pastix_complex64_t *b = (pastix_complex64_t*)(bp->b);
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
    pastix_complex32_t *sgrad = NULL;
    pastix_complex64_t  v1, v2, w;
    double normb, normx, normr, alpha, beta;
    double resid_b, eps;

    memset( &solver, 0, sizeof(struct z_solver) );
    z_refine_init( &solver, pastix_data );

    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        precond = 0;
    }

    n       = pastix_data->bcsc->n;
    itermax = pastix_data->iparm[IPARM_ITERMAX];
    eps     = pastix_data->dparm[DPARM_EPSILON_REFINEMENT];

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

    /* Allocating a vector at half-precision, NULL pointer otherwise */
    if ( pastix_data->iparm[IPARM_MIXED] )
    {
        sgrad = solver.malloc( n * sizeof(pastix_complex32_t) );
    }

    clockInit(refine_clk);clockStart(refine_clk);

    normb = solver.norm( pastix_data, n, b );
    if ( normb == 0. ) {
        normb = 1;
    }
    normx = solver.norm( pastix_data, n, x );

    /* r = b - Ax */
    solver.copy( pastix_data, n, b, gradr );
    if ( normx > 0. ) {
        solver.spmv( pastix_data, PastixNoTrans, -1., x, 1., gradr );
    }
    normr = solver.norm( pastix_data, n, gradr );

    /* r2 = r */
    solver.copy( pastix_data, n, gradr, gradr2 );
    /* p = r */
    solver.copy( pastix_data, n, gradr, gradp );

    /* resid_b = ||r|| / ||b|| */
    resid_b = normr / normb;

    while ((resid_b > eps) && (nb_iter < itermax))
    {
        clockStop((refine_clk));
        t0 = clockGet();
        nb_iter++;

        /* y = M-1 * p */
        solver.copy( pastix_data, n, gradp, grady );
        if ( precond ) {
            solver.spsv( pastix_data, grady, sgrad );
        }

        /* v = Ay */
        solver.spmv( pastix_data, PastixNoTrans, 1.0, grady, 0., gradv );

        /* alpha = (r, r2) / (v, r2) */
        alpha = solver.dot( pastix_data, n, gradv, gradr2 );
        beta  = solver.dot( pastix_data, n, gradr, gradr2 );
        alpha = beta / alpha;

        /* s = r - alpha * v */
        solver.copy( pastix_data, n, gradr, grads );
        solver.axpy( pastix_data, n, -alpha, gradv, grads );

        /* z = M^{-1} s */
        solver.copy( pastix_data, n, grads, gradz );
        if ( precond ) {
            solver.spsv( pastix_data, gradz, sgrad );
        }

        /* t = Az */
        solver.spmv( pastix_data, PastixNoTrans, 1.0, gradz, 0., gradt );

        /* w = (M-1t, M-1s) / (M-1t, M-1t) */
        /* grad2 = M-1t */
        solver.copy( pastix_data, n, gradt, grad2 );
        if ( precond ) {
            solver.spsv( pastix_data, grad2, sgrad );
        }

        /* v1 = (M-1t, M-1s) */
        /* v2 = (M-1t, M-1t) */
        v1 = solver.dot( pastix_data, n, grad2, gradz );
        v2 = solver.dot( pastix_data, n, grad2, grad2 );
        w = v1 / v2;

        /* x = x + alpha * y + w * z */
        /* x = x + alpha * y */
        solver.axpy( pastix_data, n, alpha, grady, x );

        /* x = x + w * z */
        solver.axpy( pastix_data, n, w, gradz, x );

        /* r = s - w * t*/
        solver.copy( pastix_data, n, grads, gradr );
        solver.axpy( pastix_data, n, -w, gradt, gradr );

        /* beta = (r', r2) / (r, r2) * (alpha / w) */
        /* v1 = (r', r2) */
        v1 = solver.dot( pastix_data, n, gradr, gradr2 );
        v2 = alpha / w;

        beta = v1 / beta;
        beta = beta * v2;

        /* p = r + beta * (p - w * v) */
        /* p = p - w * v */
        solver.axpy( pastix_data, n, -w, gradv, gradp );

        /* p = r + beta * p */
        solver.scal( pastix_data, n, beta, gradp );
        solver.axpy( pastix_data, n, 1., gradr, gradp );

        normr = solver.norm( pastix_data, n, gradr );
        resid_b = normr / normb;

        clockStop((refine_clk));
        t3 = clockGet();
        if ( ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) &&
             ( pastix_data->procnum == 0 ) ) {
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
    solver.free((void*) sgrad);

    return nb_iter;
}
