/**
 * @file starpu_zkernels.h
 * @precisions normal z -> s d c
 */

void starpu_zgetrfsp1d_getrf_cpu(void * buffers[], void * _args);
void starpu_zgetrfsp1d_trsm_cpu(void * buffers[], void * _args);
void starpu_zgetrfsp1d_cpu(void * buffers[], void * _args);
void starpu_zgetrfsp1d_gemm_cpu(void * buffers[], void * _args);

void starpu_zpotrfsp1d_potrf_cpu(void * buffers[], void * _args);
void starpu_zpotrfsp1d_trsm_cpu(void * buffers[], void * _args);
void starpu_zpotrfsp1d_cpu(void * buffers[], void * _args);
void starpu_zpotrfsp1d_gemm_cpu(void * buffers[], void * _args);

void starpu_zsytrfsp1d_sytrf_cpu(void * buffers[], void * _args);
void starpu_zsytrfsp1d_trsm_cpu(void * buffers[], void * _args);
void starpu_zsytrfsp1d_cpu(void * buffers[], void * _args);
void starpu_zsytrfsp1d_gemm_cpu(void * buffers[], void * _args);

void starpu_zhetrfsp1d_hetrf_cpu(void * buffers[], void * _args);
void starpu_zhetrfsp1d_trsm_cpu(void * buffers[], void * _args);
void starpu_zhetrfsp1d_cpu(void * buffers[], void * _args);
void starpu_zhetrfsp1d_gemm_cpu(void * buffers[], void * _args);
