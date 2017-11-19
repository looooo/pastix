/**
 *
 * @file pastix_parsec_gpu.c
 *
 * @author Mathieu Faverge
 * @date 2011-10-17
 * @precisions normal z -> c d s
 *
 **/
#define _GNU_SOURCE
#include "common.h"
#include "flops.h"
#include "blend/solver.h"
#include "sopalin/parsec/pastix_parsec.h"
#include "sopalin/parsec/pastix_parsec_gpu.h"

#include <parsec.h>
#include <parsec/sys/atomic.h>
#include <parsec/devices/device.h>
#include <parsec/devices/cuda/dev_cuda.h>
#include <parsec/data_internal.h>

#define PARSEC_MAX_DEVICES 32
volatile int      parsec_nbtasks_on_gpu[PARSEC_MAX_DEVICES] = { 0 };
volatile uint64_t pastix_nbflops_on_dev[PARSEC_MAX_DEVICES] = { 0 };

/* Bandwith of PCI express 16x in Go/s */
double bandwidth = 1. / 15.754e9;

double coefgemm2[][8] = {
    /* CPU Westmere Intel Xeon X5650 1D */
    { 1.94580466e-01, 3.39405974e+00, 1.50952242e+00, 3.32231901e-01,
      0.00000000e+00, 7.37678269e+00, 2.30900584e+02, 0.00000000e+00 },
    { 4.52259948e-03, 0.00000000e+00, 5.14113557e-01, 1.64044887e+00,
      9.59880994e+02, 1.40656537e+02, 8.36227071e+02, 0.00000000e+00 },
    /* CPU Westmere Intel Xeon X5650 2D */
    { 1.94580466e-01, 3.39405974e+00, 1.50952242e+00, 3.32231901e-01,
      0.00000000e+00, 7.37678269e+00, 2.30900584e+02, 0.00000000e+00 },
    { 4.52259948e-03, 0.00000000e+00, 5.14113557e-01, 1.64044887e+00,
      9.59880994e+02, 1.40656537e+02, 8.36227071e+02, 0.00000000e+00 },
};

double coefgemm[][8] =
{
    //LM Model
    //M*N*K         M*N            K*M            K*N              K              M              N             1
    /* { 4.113131e-02, 7.774065e+00,  4.266597e+00,  1.329345e+01, -7.537044e+02, -1.123117e+02, -1.503259e+03, 6.252870e+04 }, */
    /* { 7.970389e-03, 1.188658e+00,  8.456961e-02,  6.198189e-01,  2.061833e+02,  5.200194e+01,  2.901745e+02, 2.106563e+04 }, */
    //NNLS Model
    //M*N*K         M*N            K*M            K*N              K              M              N             1
    { 0.04553508,   5.52053895,    3.56445373,    6.82230469,    0.,            0.,            0.,           0.           },
    { 7.970389e-03, 1.188658e+00,  8.456961e-02,  6.198189e-01,  2.061833e+02,  5.200194e+01,  2.901745e+02, 2.106563e+04 },
};
/* { */
/*     {0.04553508,      6.82230469,     3.56445373,     5.52053895,     0.,             0.,             0.,             0.}, */
/*     {7.97038909e-03,  1.18865771e+00, 8.45696103e-02, 6.19818874e-01, 5.20019440e+01, 2.90174451e+02, 2.06183314e+02, 2.10656286e+04}, */
/* }; */

/**
 *
 * TODO: Let's not forget to add some documentation in the final version
 */
int
pastix_parsec_selectgpu_gemm1d_3p_fct( void    *arg,
                                       double   ratio,
                                       double   cpu_cost,
                                       double   gpu_cost,
                                       pastix_int_t brownum )
{
    typedef struct local_parsec_data_s {
        parsec_data_pair_t _f_C;
        parsec_data_pair_t _f_A;
        parsec_data_pair_t _f_B;
        parsec_data_pair_t unused[MAX_LOCAL_COUNT-3];
    } local_parsec_data_t;

    typedef struct local_parsec_task_s {
        PARSEC_MINIMAL_EXECUTION_CONTEXT
#if defined(PARSEC_PROF_TRACE)
        parsec_profile_data_collection_info_t prof_info;
#endif /* defined(PARSEC_PROF_TRACE) */
        assignment_t               locals[MAX_LOCAL_COUNT];
#if defined(PINS_ENABLE)
        int                        creator_core;
        int                        victim_core;
#endif /* defined(PINS_ENABLE) */
#if defined(PARSEC_SIM)
        int                        sim_exec_date;
#endif
        local_parsec_data_t         data;
    } local_parsec_task_t;

    local_parsec_task_t *this_task = (local_parsec_task_t*)arg;
    pastix_int_t dev_index, gpuid = -2;

    if (brownum == -1) {
      return -2;
    }

    dev_index = this_task->data._f_C.data_in->original->owner_device;
    if (dev_index > 1) {
        gpuid = dev_index - 2;
    }

    /**
     * If not already on the GPU, find the best one
     */
    if ( gpuid == -2 ) {
        double mintime;
        int dev, nbdevices;
        double cost;
        parsec_data_t *data;

        mintime  = cpu_cost;
        mintime  = mintime  > 0. ? mintime  : 0.;
        gpu_cost = gpu_cost > 0. ? gpu_cost : 0.;

        nbdevices = parsec_devices_enabled();
        for( dev = 2; dev < nbdevices; dev++ ) {
            cost = ( 1 + parsec_nbtasks_on_gpu[dev-2] ) * gpu_cost;

            /* Compute the transfer cost of A */
            data = this_task->data._f_A.data_in->original;
            if ( (ratio < 5) && (data->device_copies[dev] == NULL) ) {
                cost += bandwidth * data->nb_elts;
            }

            /* Compute the transfer cost of B */
            data = this_task->data._f_B.data_in->original;
            if ( (ratio < 5) && (data->device_copies[dev] == NULL) ) {
                cost += bandwidth * data->nb_elts;
            }

            /* Compute the transfer cost of B */
            data = this_task->data._f_C.data_in->original;
            if ( data->device_copies[dev] == NULL ) {
                cost += bandwidth * data->nb_elts;
            }

	    if ( cost < mintime )
            {
		gpuid = dev-2;
		mintime = cost;
            }
        }
    }

    if ( gpuid >= 0 ) {
        pastix_atomic_inc_32b( &parsec_nbtasks_on_gpu[gpuid] );
    }

    return gpuid;
}
