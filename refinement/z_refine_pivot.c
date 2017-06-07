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
 * z_grad_smp - Refine the solution using static pivoting method.
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
void z_pivot_smp (pastix_data_t *pastix_data, void *x, void *b)
{
    struct z_solver solveur;
    memset( &solveur, 0, sizeof(struct z_solver) );

    z_Pastix_Solveur(&solveur);

    pastix_bcsc_t      * bcsc           = pastix_data->bcsc;
    pastix_int_t         n              = bcsc->gN;
    Clock                refine_clk;
    double               t0             = 0;
    double               t1             = 0;
    double               t2             = 0;
    double               t3             = 0;
    pastix_complex64_t * volatile  lub  = NULL;
    pastix_complex64_t * volatile  lur  = NULL;
    pastix_complex64_t * volatile  lur2 = NULL;
    double               tmp_berr       = 0.0;
    double               berr           = 0.0;
    double               lberr          = 0.0;
    double               rberror        = 0.0;
    int                  iter           = 0;
    int                  flag           = 1;
    pastix_int_t         refinenbr        = 0.0;
    pastix_int_t         itermax;
    double               epsilonrefine;

    (void) rberror;
    (void) t0;
    (void) t1;
    (void) t2;
    itermax     = solveur.Itermax(pastix_data);
    epsilonrefine = solveur.Eps(pastix_data);

    if (pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot)
    {
        fprintf(stdout, OUT_ITERREFINE_PIVOT);
    }
    lub  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    lur  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    lur2 = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));

    solveur.B(b, lub, n);

    clockInit(refine_clk);clockStart(refine_clk);

    while(flag)
    {
        iter++;
        clockStop((refine_clk));
        t0 = clockGet();

        /* r=b-ax */
        solveur.bMAx(bcsc,lub,x,lur);
        //       z_CscbMAx(sopalin_data, me, lur, lub, sopalin_data->sopar->cscmtx,
        //               &(datacode->updovct), datacode, pastix_comm,
        //               sopar->iparm[IPARM_TRANSPOSE_SOLVE]);


        /* r'=|A||x|+|b| */
        z_bcscAxpb(PastixNoTrans, bcsc, (void *)x, (void *)lub, (void *)lur2);
        //       z_CscAxPb( sopalin_data, me, lur2, lub, sopalin_data->sopar->cscmtx,
        //                &(datacode->updovct), datacode, pastix_comm,
        //                sopar->iparm[IPARM_TRANSPOSE_SOLVE]);



        /* tmp_berr =  max_i(|lur_i|/|lur2_i|)*/
        tmp_berr = z_bcscBerr((void *)lur,(void *)lur2,n);
        //       z_CscBerr(sopalin_data, me, lur, lur2, UPDOWN_SM2XSZE,
        //               1, &tmp_berr , pastix_comm);

        berr = tmp_berr;
        if (lberr == 0)
            /* force le premier refineinement */
            lberr = 3*berr;

        print_debug(DBG_REFINE_PIVOT, "REFINE : berr lberr %6g %6g\n",
                    berr, lberr);

        /* Calcul de ||r|| et ||r||/||b|| */
        tmp_berr = z_bcscNormErr((void *)lur, (void *)lub, n);

        rberror = tmp_berr;
        print_debug(DBG_REFINE_PIVOT, "REFINE : rberror %6g\n", rberror);

        if ((refinenbr < itermax)
            && (berr > epsilonrefine)
            && (berr <= (lberr/2)))
        {

            /* LU dx = r */
            /* lur2 <= updo_vect (ie X_i)
             * updo_vect <= lur (ie B-AX_i)
             */
            memcpy(lur2, x, n * sizeof( pastix_complex64_t ));
            memcpy(x, lur, n * sizeof( pastix_complex64_t ));

            clockStop((refine_clk));
            t1 = clockGet();

            //           z_up_down_smp(arg);
            solveur.Precond(pastix_data, b, x);

            clockStop((refine_clk));
            t2 = clockGet();

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

        clockStop((refine_clk));
        t3 = clockGet();

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
        t0 = t3;
    }

    memFree_null(lub);
    memFree_null(lur);
    memFree_null(lur2);
    itermax = refinenbr;

    //   if (sopar->iparm[IPARM_END_TASK] >= PastixTaskRefine)
    //     {
    //         MUTEX_LOCK(&(sopalin_data->mutex_comm));
    //         sopalin_data->step_comm = COMMSTEP_END;
    //         print_debug(DBG_THCOMM, "%s:%d END\n", __FILE__, __LINE__);
    //         MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
    //         pthread_cond_broadcast(&(sopalin_data->cond_comm));
    //     }

    clockStop((refine_clk));
    print_debug(DBG_SOPALIN_REFINE, "%d : refinement time %lf\n", (int)me, clockGet());
    pastix_data->dparm[DPARM_REFINE_TIME] = clockGet();
}
