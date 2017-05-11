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

#if defined(PASTIX_WITH_CUDA)
#include "parsec/devices/device.h"
#include "parsec/devices/cuda/dev_cuda.h"
#endif

volatile int parsec_nbtasks_on_gpu[32] = {0};
volatile uint64_t pastix_nbflops_on_device[32] = {0};

/* Bandwith of PCI express 16x in Go/s */
double bandwidth = 15.754;

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
