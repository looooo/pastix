/**
 *
 * @file z_refine_pivot.c
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
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "bcsc.h"
#include "z_bcsc.h"
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
 * @param[out] x
 *          The solution vector.
 *
 * @param[in] b
 *          The right hand side member (only one).
 *
 *******************************************************************************
 *
 * @return Number of iterations
 *
 *******************************************************************************/
pastix_int_t z_pivot_smp (pastix_data_t *pastix_data, void *x, void *b)
{
    struct z_solver     solver;
    pastix_int_t        n;
    Clock               refine_clk;
    int                 itermax;
    pastix_int_t        iter      = 0;
    int                 precond   = 1;
    pastix_complex64_t *lur;
    pastix_complex64_t *lur2;
    double eps;

    pastix_bcsc_t       *bcsc = pastix_data->bcsc;
    double               tmp_berr       = 0.0;
    double               berr           = 0.0;
    double               lberr          = 0.0;
    int                  flag           = 1;
    pastix_int_t         refinenbr     = 0.0;

    memset( &solver, 0, sizeof(struct z_solver) );
    z_Pastix_Solver(&solver);

    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        precond = 0;
    }

    n       = solver.getN(    pastix_data );
    itermax = solver.getImax( pastix_data );
    eps     = solver.getEps(  pastix_data );

    if (pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot)
    {
        fprintf(stdout, OUT_ITERREFINE_PIVOT);
    }
    lur  = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    lur2 = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));

    clockInit(refine_clk);clockStart(refine_clk);

    while(flag)
    {
        iter++;
        clockStop((refine_clk));

        /* Compute r = b - A * x */
        solver.copy( n, b, lur );
        solver.spmv( pastix_data, -1., x, 1., lur );

        /* r' = |A||x| + |b| */
        z_bcscAxpb(PastixNoTrans, bcsc, x, b, lur2);

        /* tmp_berr =  max_i(|lur_i|/|lur2_i|)*/
        tmp_berr = z_bcscBerr( lur, lur2, n );

        berr = tmp_berr;
        if (lberr == 0.) {
            /* Force te first error */
            lberr = 3*berr;
        }

        /* Compute ||r|| and ||r||/||b|| */
        tmp_berr = z_bcscNormErr((void *)lur, (void *)b, n);

        if ((refinenbr < itermax)
            && (berr > eps)
            && (berr <= (lberr/2)))
        {

            /* LU dx = r */
            /* lur2 <= updo_vect (ie X_i)
             * updo_vect <= lur (ie B-AX_i)
             */
            memcpy(lur2, x, n * sizeof( pastix_complex64_t ));
            memcpy(x, lur, n * sizeof( pastix_complex64_t ));

            clockStop((refine_clk));

            //           z_up_down_smp(arg);
            solver.copy( n, b, x );
            if ( precond ) {
                solver.spsv( pastix_data, x );
            }

            clockStop(refine_clk);

            /* updo_vect <= updo_vect (ie PRECOND(B-AX_i)) + lur2 (ie X_i) */
            z_bcscAxpy( n, 1, 1.0, (void*)lur2, x );

            /* lastberr = berr */
            lberr = berr;
            refinenbr++;
        }
        else
        {
            flag = 0;
        }

        clockStop(refine_clk);

        //       if (sopar->iparm[IPARM_VERBOSE] > PastixVerboseNot)
        //         {
        //           double sst, rst = 0.0;
        //           double stt, rtt;
        //           double err, berr = sopalin_data->berr;
        //
        //           stt = t3 - t0;
        //           sst = t2-t1;
        //           MyMPI_Reduce(&sst, &rst, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
        //
        //           MyMPI_Reduce(&berr, &err, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
        //           MyMPI_Reduce(&stt,  &rtt, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
        //           if (SOLV_PROCNUM == 0)
        //             {
        //               fprintf(stdout, OUT_ITERREFINE_ITER, (int)sopalin_data->refinenbr);
        //               fprintf(stdout, OUT_ITERREFINE_TTS, rst);
        //               fprintf(stdout, OUT_ITERREFINE_TTT, rtt);
        //               fprintf(stdout, OUT_ITERREFINE_ERR, err);
        //             }
        //         }
    }

    memFree_null(lur);
    memFree_null(lur2);
    itermax = refinenbr;

    //   if (sopar->iparm[IPARM_END_TASK] >= PastixTaskRefine)
    //     {
    //         MUTEX_LOCK(&(sopalin_data->mutex_comm));
    //         sopalin_data->step_comm = COMMSTEP_END;
    //         MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
    //         pthread_cond_broadcast(&(sopalin_data->cond_comm));
    //     }

    clockStop((refine_clk));
    pastix_data->dparm[DPARM_REFINE_TIME] = clockGet();

    return iter;
}
