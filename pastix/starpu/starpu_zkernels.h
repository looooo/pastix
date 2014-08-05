/**
 * @file starpu_zkernels.h
 * @precisions normal z -> s d c
 */
#ifndef STARPU_ZKERNELS_H
#define STARPU_ZKERNELS_H

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

void starpu_zgetrfsp1d_getrf_cpu(void * buffers[], void * _args);
void starpu_zgetrfsp1d_trsm_cpu(void * buffers[], void * _args);
void starpu_zgetrfsp1d_cpu(void * buffers[], void * _args);
void starpu_zgetrfsp1d_gemm_cpu(void * buffers[], void * _args);
void starpu_zgetrfsp1d_geadd_cpu(void * buffers[], void * _args);

void starpu_zpotrfsp1d_potrf_cpu(void * buffers[], void * _args);
void starpu_zpotrfsp1d_trsm_cpu(void * buffers[], void * _args);
void starpu_zpotrfsp1d_cpu(void * buffers[], void * _args);
void starpu_zpotrfsp1d_gemm_cpu(void * buffers[], void * _args);
void starpu_zpotrfsp1d_poadd_cpu(void * buffers[], void * _args);

void starpu_zsytrfsp1d_sytrf_cpu(void * buffers[], void * _args);
void starpu_zsytrfsp1d_trsm_cpu(void * buffers[], void * _args);
void starpu_zsytrfsp1d_cpu(void * buffers[], void * _args);
void starpu_zsytrfsp1d_gemm_cpu(void * buffers[], void * _args);
void starpu_zsytrfsp1d_gemm_cuda(void * buffers[], void * _args);
void starpu_zsytrfsp1d_syadd_cpu(void * buffers[], void * _args);

void starpu_zhetrfsp1d_hetrf_cpu(void * buffers[], void * _args);
void starpu_zhetrfsp1d_trsm_cpu(void * buffers[], void * _args);
void starpu_zhetrfsp1d_cpu(void * buffers[], void * _args);
void starpu_zhetrfsp1d_gemm_cpu(void * buffers[], void * _args);
void starpu_zhetrfsp1d_headd_cpu(void * buffers[], void * _args);

#  ifndef STARPU_PROF_CALLBACK_DEFINE
#  define STARPU_PROF_CALLBACK_DEFINE
void starpu_prof_callback(void *callback_arg);
#  endif

/* For old kernels ===> TRASH when not used anymore*/
#define ARCH_CPU  0
#define ARCH_CUDA 1


#define trsm_starpu_common                API_CALL(z_trsm_starpu_common)
#define xxtrf_starpu_common               API_CALL(z_xxtrf_starpu_common)
#define trfsp1d_starpu_common             API_CALL(z_trfsp1d_starpu_common)
#define trfsp1d_gemm_starpu_common        API_CALL(z_trfsp1d_gemm_starpu_common)
#define trfsp1d_sparse_gemm_starpu_common API_CALL(z_trfsp1d_sparse_gemm_starpu_common)
#define trsm_starpu_cpu                   API_CALL(z_trsm_starpu_cpu)
#define xxtrf_starpu_cpu                  API_CALL(z_xxtrf_starpu_cpu)
#define trfsp1d_starpu_cpu                API_CALL(z_trfsp1d_starpu_cpu)
#define trfsp1d_gemm_starpu_cpu           API_CALL(z_trfsp1d_gemm_starpu_cpu)
#define trfsp1d_sparse_gemm_starpu_cpu    API_CALL(z_trfsp1d_sparse_gemm_starpu_cpu)
#define trsm_starpu_cuda                  API_CALL(z_trsm_starpu_cuda)
#define xxtrf_starpu_cuda                 API_CALL(z_xxtrf_starpu_cuda)
#define trfsp1d_starpu_cuda               API_CALL(z_trfsp1d_starpu_cuda)
#define trfsp1d_gemm_starpu_cuda          API_CALL(z_trfsp1d_gemm_starpu_cuda)
#define trfsp1d_sparse_gemm_starpu_cuda   API_CALL(z_trfsp1d_sparse_gemm_starpu_cuda)

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

#endif
