/**
 *
 * @file z_refine_pivot.c
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Theophile Terraz
 * @author Vincent Bridonneau
 * @date 2024-07-05
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "z_refine_functions.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_refine
 *
 * z_pivot_smp - Refine the solution using static pivoting method.
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
z_pivot_smp( pastix_data_t *pastix_data,
             pastix_rhs_t   xp,
             pastix_rhs_t   bp )
{
    struct z_solver     solver;
    Clock               t0, t3, refine_clk;
    pastix_int_t        n, itermax;
    pastix_complex64_t *x = (pastix_complex64_t*)(xp->b);
    pastix_complex64_t *b = (pastix_complex64_t*)(bp->b);
    pastix_complex64_t *r, *dx;
    pastix_complex32_t *sb = NULL;
    double              eps, normb, normr;
    double              berr, last_berr;
    pastix_int_t        iter = 0;
    int                 flag = 1;

    memset( &solver, 0, sizeof(struct z_solver) );
    z_refine_init( &solver, pastix_data );

    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        fprintf(stderr, "pastix_task_refine: Simple refinement cannot be applied without preconditionner\n" );
        return -1;
    }

    n       = pastix_data->bcsc->n;
    itermax = pastix_data->iparm[IPARM_ITERMAX];
    eps     = pastix_data->dparm[DPARM_EPSILON_REFINEMENT];

    if (pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot)
    {
        fprintf(stdout, OUT_ITERREFINE_PIVOT);
    }
    r  = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    dx = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));

    clockInit(refine_clk);
    clockStart(refine_clk);

    normb = solver.norm( pastix_data, n, b );
    if ( normb == 0. ) {
        normb = 1;
    }

    /* Allocating a vector at half-precision, NULL pointer otherwise */
    if ( pastix_data->iparm[IPARM_MIXED] )
    {
        sb = solver.malloc( n * sizeof(pastix_complex32_t) );
    }

    t0 = clockGet();
    while(flag)
    {
        /* Compute r = b - A * x */
        solver.copy( pastix_data, n, b, r );
        solver.spmv( pastix_data, PastixNoTrans, -1., x, 1., r );

        /*
         * berr should be equal to the componentwise backward error in the literature:
         *     max( r_i / ( |A| |x| + |b| )_i )
         * For simplicity, we replace it by ||r||_f / ||b||_f which may not be
         * as good as the previous one.
         */
        normr = solver.norm( pastix_data, n, r );
        berr = normr / normb;

        /* Force the first error */
        if ( iter == 0 ) {
            last_berr = 3 * berr;
        }
        else {
            t3 = clockGet();
            if ( ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) &&
                 ( pastix_data->procnum == 0 ) ) {
                solver.output_oneiter( t0, t3, berr, iter );
            }
            t0 = clockGet();
        }

        if ( (iter < itermax) &&
             (berr > eps) &&
             (berr <= (last_berr / 2.)) )
        {
            t3 = clockGet();

            /* Solve A dx = r */
            solver.copy( pastix_data, n, r, dx );
            solver.spsv( pastix_data, dx, sb );

            /* Accumulate the solution: x = x + dx */
            solver.axpy( pastix_data, n, 1.0, dx, x );

            last_berr = berr;
        }
        else
        {
            flag = 0;
        }

        iter++;
    }
    clockStop(refine_clk);

    solver.output_final( pastix_data, berr, iter, refine_clk, x, x );

    solver.free( r );
    solver.free( dx );
    solver.free( sb );

    return iter;
}
