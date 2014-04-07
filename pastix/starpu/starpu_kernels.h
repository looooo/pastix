#ifndef STARPU_KERNELS_H
#define STARPU_KERNELS_H

#include "redefine_functions.h"
#include "sopalin_define.h"

#define ARCH_CPU  0
#define ARCH_CUDA 1

#define MAT_zaxpy( m, n, alpha,                         \
                   A, lda,                              \
                   B, ldb )                             \
    {                                                   \
        pastix_int_t i,j;                                 \
        for (j=0; j<n; j++) {                           \
            for (i=0; i<m; i++) {                       \
                B[j*ldb+i] = B[j*ldb+i] - A[j*lda+i];   \
            }                                           \
        }                                               \
    }

#define MAT_zaxpyt( m, n, alpha,                        \
                    A, lda,                             \
                    B, ldb )                            \
    {                                                   \
        pastix_int_t i,j;                                 \
        for (j=0; j<n; j++) {                           \
            for (i=0; i<m; i++) {                       \
                B[j*ldb+i] = B[j*ldb+i] - A[j+i*lda];   \
            }                                           \
        }                                               \
    }


#define trsm_starpu_common                API_CALL(trsm_starpu_common)
#define xxtrf_starpu_common               API_CALL(xxtrf_starpu_common)
#define trfsp1d_starpu_common             API_CALL(trfsp1d_starpu_common)
#define trfsp1d_gemm_starpu_common        API_CALL(trfsp1d_gemm_starpu_common)
#define trfsp1d_sparse_gemm_starpu_common API_CALL(trfsp1d_sparse_gemm_starpu_common)
#define trsm_starpu_cpu                   API_CALL(trsm_starpu_cpu)
#define xxtrf_starpu_cpu                  API_CALL(xxtrf_starpu_cpu)
#define trfsp1d_starpu_cpu                API_CALL(trfsp1d_starpu_cpu)
#define trfsp1d_gemm_starpu_cpu           API_CALL(trfsp1d_gemm_starpu_cpu)
#define trfsp1d_sparse_gemm_starpu_cpu    API_CALL(trfsp1d_sparse_gemm_starpu_cpu)
#define trsm_starpu_cuda                  API_CALL(trsm_starpu_cuda)
#define xxtrf_starpu_cuda                 API_CALL(xxtrf_starpu_cuda)
#define trfsp1d_starpu_cuda               API_CALL(trfsp1d_starpu_cuda)
#define trfsp1d_gemm_starpu_cuda          API_CALL(trfsp1d_gemm_starpu_cuda)
#define trfsp1d_sparse_gemm_starpu_cuda   API_CALL(trfsp1d_sparse_gemm_starpu_cuda)

void trfsp1d_starpu_cpu(void * buffers[], void * _args);
void xxtrf_starpu_cpu(void * buffers[], void * _args);
void trsm_starpu_cpu(void * buffers[], void * _args);
void trfsp1d_gemm_starpu_cpu(void * buffers[], void * _args);
void trfsp1d_sparse_gemm_starpu_cpu(void * buffers[], void * _args);
void trfsp1d_starpu_cuda(void * buffers[], void * _args);
void xxtrf_starpu_cuda(void * buffers[], void * _args);
void trsm_starpu_cuda(void * buffers[], void * _args);
void trfsp1d_gemm_starpu_cuda(void * buffers[], void * _args);
void trfsp1d_sparse_gemm_starpu_cuda(void * buffers[], void * _args);

#define init_coeftab_cpu API_CALL(init_coeftab_cpu)
#define init_coeftab_cuda API_CALL(init_coeftab_cuda)
void init_coeftab_cpu(void * buffers[], void * _args);
void init_coeftab_cuda(void * buffers[], void * _args);

#define fill_coeftab_cpu API_CALL(fill_coeftab_cpu)
#define fill_coeftab_cuda API_CALL(fill_coeftab_cuda)
void fill_coeftab_cpu(void * buffers[], void * _args);
void fill_coeftab_cuda(void * buffers[], void * _args);
#endif /* STARPU_KERNELS_H */
