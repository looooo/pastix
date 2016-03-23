/**
 *
 * @file core_zgelrops.c
 *
 *  PaStiX kernel routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Grégoire Pichon
 * @date 2016-23-03
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "pastix_zcores.h"
#include <cblas.h>
#include <lapacke.h>

static pastix_complex64_t zone  =  1.;
static pastix_complex64_t mzone = -1.;
static pastix_complex64_t zzero =  0.;

pastix_int_t z_compress_LR(pastix_complex64_t *fL,
                           pastix_int_t stride,
                           pastix_int_t dimb,
                           pastix_int_t dima,
                           pastix_complex64_t *u,
                           pastix_int_t ldu,
                           pastix_complex64_t *v,
                           pastix_int_t ldv){

    pastix_complex64_t *block;
    pastix_int_t        ret;
    pastix_int_t        i;
    double             *s;
    double             *superb;

    pastix_int_t dim_min = dima;
    if (dimb < dim_min)
        dim_min = dimb;

    /* Note that we have to copy fL because LAPACKE_zgesvd erases the matrix */
    block  = malloc( dima * dimb * sizeof(pastix_complex64_t));
    s      = malloc( dim_min * sizeof(double));
    superb = malloc( dimb * sizeof(double));

    for (i=0; i<dima; i++){
        memcpy( block + i * dimb, fL + i * stride, dimb * sizeof(pastix_complex64_t));
    }

    ret = LAPACKE_zgesvd( CblasColMajor, 'S', 'S',
                          dimb, dima, block, dimb,
                          s, u, ldu, v, ldv, superb );

    if( ret != 0 ){
        printf("SVD FAILED %ld\n\n", ret);
        exit(1);
    }

    char *tol        = getenv("TOLERANCE");
    double tolerance = atof(tol);

    pastix_int_t rank = dim_min;
    for (i=0; i<dim_min-1; i++){
        if (s[i] / s[0] < tolerance){
            rank = i+1;
            break;
        }
    }

    for (i=0; i<rank; i++){
        cblas_dscal(dimb, s[i], &(u[i*ldu]), 1);
    }

    free(block);
    free(s);
    free(superb);
    return rank;
}
