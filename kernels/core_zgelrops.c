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
static pastix_complex64_t zzero =  0.;

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_z_compress_LR - Compresses a dense block into a u v^T LR structure.
 *
 *******************************************************************************
 *
 * @param[in] fL
 *          Pointer to the dense structure of size dimb * dima
 *          Leading dimension is stride
 *
 * @param[out] u
 *          Pointer to the u factor of LR representation of size dimb * rank
 *          Leading dimension is ldu
 *
 * @param[out] v
 *          Pointer to the v factor of LR representation of size dima * rank
 *          Leading dimension is ldv
 *          Note that due to LAPACKE_zgesvd this block is stored transposed
 *
 *
 *******************************************************************************
 *
 * @return
 *          The rank of the compressed structure.
 *
 *******************************************************************************/
pastix_int_t core_z_compress_LR(pastix_complex64_t *fL,
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
    block  = malloc( dima * dimb * sizeof(pastix_complex64_t) );
    s      = malloc( dim_min * sizeof(double) );
    superb = malloc( dimb * sizeof(double) );

    for (i=0; i<dima; i++){
        memcpy( block + i * dimb, fL + i * stride, dimb * sizeof(pastix_complex64_t) );
    }

    ret = LAPACKE_zgesvd( CblasColMajor, 'S', 'S',
                          dimb, dima,
                          block, dimb,
                          s, u, ldu, v, ldv,
                          superb );

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

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_z_uncompress_LR - Uncompresses a u v^T LR structure into a dense block.
 *
 *******************************************************************************
 *
 * @param[out] fL
 *          Pointer to the dense structure of size dimb * dima
 *          Leading dimension is stride
 *
 * @param[in] u
 *          Pointer to the u factor of LR representation of size dimb * rank
 *          Leading dimension is ldu
 *
 * @param[in] v
 *          Pointer to the v factor of LR representation of size dima * rank
 *          Leading dimension is ldv
 *          Note that due to LAPACKE_zgesvd this block is stored transposed
 *
 *
 *******************************************************************************
 *
 * @return
 *          The rank of the compressed structure.
 *
 *******************************************************************************/
void core_z_uncompress_LR(pastix_complex64_t *fL,
                          pastix_int_t stride,
                          pastix_int_t dimb,
                          pastix_int_t dima,
                          pastix_complex64_t *u,
                          pastix_int_t ldu,
                          pastix_complex64_t *v,
                          pastix_int_t ldv,
                          pastix_int_t rank){

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                dimb, dima, rank,
                CBLAS_SADDR(zone),  u,  ldu,
                                    v,  ldv,
                CBLAS_SADDR(zzero), fL, stride);
}
