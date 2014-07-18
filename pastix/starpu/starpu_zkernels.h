/**
 * @file starpu_zkernels.h
 * @precisions normal z -> s d c
 */
#ifndef STARPU_ZKERNELS_H
#define STARPU_ZKERNELS_H

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
#endif
