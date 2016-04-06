/**
 *
 * @file test_zgetrfsp.c
 *
 *  Testing for PaStiX kernel routines
 *  PaStiX is a asoftware package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/

#include <malloc.h>
#include <float.h>
#include "common.h"
#include "core_zgetrfsp.c" //#include "pastix_zcores.h"

int
test_zgetf2sp(pastix_int_t m, pastix_int_t n, pastix_int_t lda) {
    pastix_int_t nbpivot = 0;
    pastix_int_t criteria = 1e-12;
    pastix_int_t i,j,k;

    pastix_complex64_t *A    = malloc(n*lda*sizeof(pastix_complex64_t));
    pastix_complex64_t *LU   = malloc(n*lda*sizeof(pastix_complex64_t));
    pastix_complex64_t *Acpy = malloc(n*lda*sizeof(pastix_complex64_t));
    double diff, norm, result;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            double a=((double)rand()/(double)RAND_MAX);
            double b=((double)rand()/(double)RAND_MAX);
            A[i*lda + j]    = a /*+  I*b */;
            Acpy[i*lda + j] = a /*+  I*b */;
        }
    }
    core_zgetf2sp(m, n, A, lda,
                   &nbpivot,
                   criteria );

    memset(LU, 0, n*lda*sizeof(pastix_complex64_t));
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            for ( k = 0; k < pastix_imin(i+1,j+1); k++) {
                if (k == j)
                    LU[i*lda+j] += A[i*lda+k];
                else
                    LU[i*lda+j] += A[k*lda+j]*A[i*lda+k];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            diff = (Acpy[i*lda+j] - LU[i*lda+j]) *
                conj(Acpy[i*lda+j]  - LU[i*lda+j]);
            norm = (Acpy[i*lda+j]) * conj(Acpy[i*lda+j]);
        }
    }
    result = sqrt(diff)/sqrt(norm)*n*DBL_EPSILON;
    fprintf(stdout, "||LU - A||/(||A||.N.eps) %.3g\n", result);
    if ( isnan(diff)
         || isinf(diff)
         || isnan(result)
         || isinf(result)
         || (result > 60.0) ) {
        fprintf(stdout,"-- Factorization is suspicious ! \n");
        return 1;
    }
    else{
        fprintf(stdout,"-- Factorization is CORRECT ! \n");
        return 0;
    }
}

int
test_zgetrfsp(pastix_int_t n, pastix_int_t lda) {
    pastix_int_t nbpivot = 0;
    pastix_int_t criteria = 1e-12;
    pastix_int_t i,j,k;

    pastix_complex64_t *A    = malloc(n*lda*sizeof(pastix_complex64_t));
    pastix_complex64_t *LU   = malloc(n*lda*sizeof(pastix_complex64_t));
    pastix_complex64_t *Acpy = malloc(n*lda*sizeof(pastix_complex64_t));
    double diff, norm, result;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double a=((double)rand()/(double)RAND_MAX);
            double b=((double)rand()/(double)RAND_MAX);
            A[i*lda + j]    = a /*+  I*b */;
            Acpy[i*lda + j] = a /*+  I*b */;
        }
    }
    core_zgetrfsp(n, A, lda,
                  &nbpivot,
                  criteria );

    memset(LU, 0, n*lda*sizeof(pastix_complex64_t));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for ( k = 0; k < pastix_imin(i+1,j+1); k++) {
                if (k == j)
                    LU[i*lda+j] += A[i*lda+k];
                else
                    LU[i*lda+j] += A[k*lda+j]*A[i*lda+k];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            diff = (Acpy[i*lda+j] - LU[i*lda+j]) *
                conj(Acpy[i*lda+j]  - LU[i*lda+j]);
            norm = (Acpy[i*lda+j]) * conj(Acpy[i*lda+j]);
        }
    }
    result = sqrt(diff)/sqrt(norm)*n*DBL_EPSILON;
    fprintf(stdout, "||LU - A||/(||A||.N.eps) %.3g\n", result);
    if ( isnan(diff)
         || isinf(diff)
         || isnan(result)
         || isinf(result)
         || (result > 60.0) ) {
        fprintf(stdout,"-- Factorization is suspicious ! \n");
        return 1;
    }
    else{
        fprintf(stdout,"-- Factorization is CORRECT ! \n");
        return 0;
    }
}

int main( int argc, char**argv) {
    pastix_int_t ret = 0;
    fprintf(stdout, "test_zgetf2sp M = N = LDA = 128\n");
    ret |= test_zgetf2sp(128, 128, 128);
    fprintf(stdout, "test_zgetf2sp M = N = 128, LDA = 256\n");
    ret |= test_zgetf2sp(128, 128, 256);
    fprintf(stdout, "test_zgetrfsp N = 256, LDA = 512\n");
    ret |= test_zgetrfsp(256, 512);
    return ret;
}

