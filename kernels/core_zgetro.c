/**
 *
 * @file core_zgetro.c
 *
 *  PaStiX kernel routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @copyright 2010-2015 Univ. of Tennessee, Univ. of California Berkeley and
 *  Univ. of Colorado Denver. All rights reserved.
 * @copyright 2015-2017
 *  Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights
 *  reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/**
 ******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zgetro transposes a  m-by-n matrix in place using an extra workspace of
 * size m-by-n.
 *
 *******************************************************************************
 *
 * @param[in] m
 *         Number of rows of A.
 *
 * @param[in] n
 *         Number of columns of A.
 *
 * @param[in] A
 *         Matrix to be transposed.
 *
 * @param[in] lda
 *         Leading dimension of matrix A.
 *
 * @param[in] B
 *         On exit B = trans(A).
 *
 * @param[in] ldb
 *         Leading dimension of matrix B.
 *
 ******************************************************************************/
void core_zgetro(int m, int n,
                 const pastix_complex64_t *A, int lda,
                 pastix_complex64_t *B, int ldb)
{
    int i, j;

    /* rectangular transposition (use workspace) */
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            B[j+i*ldb] = A[i+j*lda];
        }
    }
}

/* TODO: remove? apparently not used... */
/* void core_zaxpyt(int m, int n, pastix_complex64_t alpha, */
/*                  const pastix_complex64_t *A, int lda, */
/*                  pastix_complex64_t *B, int ldb) */
/* { */
/*     int i, j; */

/*     /\* rectangular transposition (use workspace) *\/ */
/*     for (j=0; j<n; j++) { */
/*         for (i=0; i<m; i++) { */
/*             B[j*ldb+i] = B[j*ldb+i] + alpha * A[j+i*lda]; */
/*         } */
/*     } */
/* } */
