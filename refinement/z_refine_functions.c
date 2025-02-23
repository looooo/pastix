/**
 *
 * @file z_refine_functions.c
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Theophile Terraz
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @author Vincent Bridonneau
 * @date 2024-07-05
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "cblas.h"
#include "bcsc/bcsc.h"
#include "bcsc/bvec.h"
#include "bcsc/bcsc_z.h"
#include "sopalin/sopalin_data.h"
#include "refinement/z_refine_functions.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Print statistics about one iteration
 *
 *******************************************************************************
 *
 * @param[in] t0
 *          The clock value at the beginning of the iteration
 *
 * @param[in] tf
 *          The clock value at the end of the iteration
 *
 * @param[in] err
 *          The backward error after the iteration
 *
 * @param[in] nb_iters
 *          Current number of refinement iterations
 *
 *******************************************************************************/
void z_refine_output_oneiter( pastix_fixdbl_t t0, pastix_fixdbl_t tf, double err, pastix_int_t nb_iters )
{
    pastix_fixdbl_t stt;

    stt = tf - t0;
    fprintf(stdout, OUT_ITERREFINE_ITER, (int)nb_iters);
    fprintf(stdout, OUT_ITERREFINE_TTT, stt);
    fprintf(stdout, OUT_ITERREFINE_ERR, err);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Final output
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[in] err
 *          The final backward error
 *
 * @param[in] nb_iters
 *          The final number of iterations
 *
 * @param[in] tf
 *          The final clock value
 *
 * @param[inout] x
 *          The vector that is to be overwritten by gmresx
 *
 * @param[in] gmresx
 *          The final solution
 *
 *******************************************************************************/
void z_refine_output_final( pastix_data_t      *pastix_data,
                            pastix_complex64_t  err,
                            pastix_int_t        nb_iters,
                            pastix_fixdbl_t     tf,
                            void               *x,
                            pastix_complex64_t *gmresx )
{
    (void)pastix_data;
    (void)err;
    (void)nb_iters;
    (void)tf;
    (void)x;
    (void)gmresx;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Initiate functions pointers to define basic operations
 *
 *******************************************************************************
 *
 * @param[out] solver
 *          The structure to be filled
 *
 * @param[in] pastix_data
 *          TODO
 *
 *******************************************************************************/
void z_refine_init( struct z_solver *solver,
                    pastix_data_t   *pastix_data )
{
    pastix_scheduler_t sched = pastix_data->iparm[IPARM_SCHEDULER];

    /* Allocations */
    solver->malloc  = &bvec_malloc;
    solver->free    = &bvec_free;

    /* Output */
    solver->output_oneiter = &z_refine_output_oneiter;
    solver->output_final   = &z_refine_output_final;

    /* Basic operations */
    solver->spsv = &bcsc_zspsv;
    if ( sched == PastixSchedSequential )
    {
        solver->spmv = &bcsc_zspmv;
        solver->copy = &bvec_zcopy_seq;
        solver->dot  = &bvec_zdotc_seq;
        solver->axpy = &bvec_zaxpy_seq;
        solver->scal = &bvec_zscal_seq;
        solver->norm = &bvec_znrm2_seq;
        solver->gemv = &bvec_zgemv_seq;
    } else {
        solver->spmv = &bcsc_zspmv;
        solver->copy = &bvec_zcopy_smp;
        solver->dot  = &bvec_zdotc_smp;
        solver->axpy = &bvec_zaxpy_smp;
        solver->scal = &bvec_zscal_smp;
        solver->norm = &bvec_znrm2_smp;
        solver->gemv = &bvec_zgemv_smp;
    }
}
