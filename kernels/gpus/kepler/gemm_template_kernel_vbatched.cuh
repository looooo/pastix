/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Mark Gates
       @author Azzam Haidar
       @author Ahmad Abdelfattah
*/

#ifndef GEMM_TEMPLATE_KERNEL_VBATCHED_CUH
#define GEMM_TEMPLATE_KERNEL_VBATCHED_CUH

static inline int pastix_ceildiv( int a, int b ) {
    return (a + b -1)/b;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
#include "gemm_template_device_defs.cuh"
#include "gemm_template_device.cuh"
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int CONJA, const int CONJB>
static __global__
void gemm_template_vbatched_nn_kernel(
    int* M, int* N, int* K,
    T const * const * Aarray, int* LDA,
    T const * const * Barray, int* LDB,
    T**       Carray, int* LDC,
    T alpha, T beta,
    int offsetA, int offsetB)
{
    const int batchid = blockIdx.z;
    const int my_M = M[batchid];
    const int my_N = N[batchid];
    const int my_K = K[batchid];
    if(my_M <= 0 || my_N <= 0 || my_K <= 0) return;
    if( Aarray[batchid] == NULL || Barray[batchid] == NULL || Carray[batchid] == NULL ) return;
    if( blockIdx.x >= (my_M+BLK_M-1)/BLK_M ) return;
    if( blockIdx.y >= (my_N+BLK_N-1)/BLK_N ) return;

    /*
    #ifdef TEXTURE_1D
    int matrixA_size = gridDim.z > 1 ?  Aarray[1] - Aarray[0] : 0;
    int matrixB_size = gridDim.z > 1 ?  Aarray[1] - Aarray[0] : 0;
    offsetA += batchid*matrixA_size;
    offsetB += batchid*matrixB_size;
    #endif
    */

    gemm_template_device_nn<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB, (BLK_M/DIM_X), (BLK_N/DIM_Y), CONJA, CONJB>
    ( my_M, my_N, my_K, Aarray[batchid], LDA[batchid], Barray[batchid], LDB[batchid], Carray[batchid], LDC[batchid], alpha, beta, offsetA, offsetB );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int CONJA, const int CONJB>
static __global__
void gemm_template_vbatched_nt_kernel(
    int* M, int* N, int* K,
    T const * const * Aarray, int* LDA,
    T const * const * Barray, int* LDB,
    T**       Carray, int* LDC,
    T alpha, T beta,
    int offsetA, int offsetB)
{
    const int batchid = blockIdx.z;
    const int my_M = M[batchid];
    const int my_N = N[batchid];
    const int my_K = K[batchid];
    if(my_M <= 0 || my_N <= 0 || my_K <= 0) return;
    if( Aarray[batchid] == NULL || Barray[batchid] == NULL || Carray[batchid] == NULL ) return;
    if( blockIdx.x >= (my_M+BLK_M-1)/BLK_M ) return;
    if( blockIdx.y >= (my_N+BLK_N-1)/BLK_N ) return;

    /*
    #ifdef TEXTURE_1D
    int matrixA_size = gridDim.z > 1 ?  Aarray[1] - Aarray[0] : 0;
    int matrixB_size = gridDim.z > 1 ?  Aarray[1] - Aarray[0] : 0;
    offsetA += batchid*matrixA_size;
    offsetB += batchid*matrixB_size;
    #endif
    */

    gemm_template_device_nt<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB, (BLK_M/DIM_X), (BLK_N/DIM_Y), CONJA, CONJB>
    ( my_M, my_N, my_K, Aarray[batchid], LDA[batchid], Barray[batchid], LDB[batchid], Carray[batchid], LDC[batchid], alpha, beta, offsetA, offsetB );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int CONJA, const int CONJB>
static __global__
void gemm_template_vbatched_tn_kernel(
    int* M, int* N, int* K,
    T const * const * Aarray, int* LDA,
    T const * const * Barray, int* LDB,
    T**       Carray, int* LDC,
    T alpha, T beta,
    int offsetA, int offsetB)
{
    const int batchid = blockIdx.z;
    const int my_M = M[batchid];
    const int my_N = N[batchid];
    const int my_K = K[batchid];
    if(my_M <= 0 || my_N <= 0 || my_K <= 0) return;
    if( Aarray[batchid] == NULL || Barray[batchid] == NULL || Carray[batchid] == NULL ) return;
    if( blockIdx.x >= (my_M+BLK_M-1)/BLK_M ) return;
    if( blockIdx.y >= (my_N+BLK_N-1)/BLK_N ) return;

    /*
    #ifdef TEXTURE_1D
    int matrixA_size = gridDim.z > 1 ?  Aarray[1] - Aarray[0] : 0;
    int matrixB_size = gridDim.z > 1 ?  Aarray[1] - Aarray[0] : 0;
    offsetA += batchid*matrixA_size;
    offsetB += batchid*matrixB_size;
    #endif
    */

    gemm_template_device_tn<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB, (BLK_M/DIM_X), (BLK_N/DIM_Y), CONJA, CONJB>
    ( my_M, my_N, my_K, Aarray[batchid], LDA[batchid], Barray[batchid], LDB[batchid], Carray[batchid], LDC[batchid], alpha, beta, offsetA, offsetB );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int CONJA, const int CONJB>
static __global__
void gemm_template_vbatched_tt_kernel(
    int* M, int* N, int* K,
    T const * const * Aarray, int* LDA,
    T const * const * Barray, int* LDB,
    T**       Carray, int* LDC,
    T alpha, T beta,
    int offsetA, int offsetB)
{
    const int batchid = blockIdx.z;
    const int my_M = M[batchid];
    const int my_N = N[batchid];
    const int my_K = K[batchid];
    if(my_M <= 0 || my_N <= 0 || my_K <= 0) return;
    if( Aarray[batchid] == NULL || Barray[batchid] == NULL || Carray[batchid] == NULL ) return;
    if( blockIdx.x >= (my_M+BLK_M-1)/BLK_M ) return;
    if( blockIdx.y >= (my_N+BLK_N-1)/BLK_N ) return;

    /*
    #ifdef TEXTURE_1D
    int matrixA_size = gridDim.z > 1 ?  Aarray[1] - Aarray[0] : 0;
    int matrixB_size = gridDim.z > 1 ?  Aarray[1] - Aarray[0] : 0;
    offsetA += batchid*matrixA_size;
    offsetB += batchid*matrixB_size;
    #endif
    */

    gemm_template_device_tt<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB, (BLK_M/DIM_X), (BLK_N/DIM_Y), CONJA, CONJB>
    ( my_M, my_N, my_K, Aarray[batchid], LDA[batchid], Barray[batchid], LDB[batchid], Carray[batchid], LDC[batchid], alpha, beta, offsetA, offsetB );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
// kernel wrappers
////////////////////////////////////////////////////////////////////////////////////////////////////
// NN
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K, const int dim_vec,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int CONJA, const int CONJB>
void gemm_template_vbatched_nn(
    pastix_int_t* m, pastix_int_t* n, pastix_int_t* k,
    T const * const * dA_array, pastix_int_t* ldda,
    T const * const * dB_array, pastix_int_t* lddb,
    T**       dC_array, pastix_int_t* lddc,
    T alpha, T beta,
    pastix_int_t offsetA, pastix_int_t offsetB,
    pastix_int_t batchCount, cudaStream_t stream,
    pastix_int_t max_m, pastix_int_t max_n)
{
    dim3 dimBlock(DIM_X, DIM_Y);
    dim3 dimGrid( pastix_ceildiv( max_m, BLK_M ), pastix_ceildiv( max_n, BLK_N ), batchCount );
    gemm_template_vbatched_nn_kernel<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB, CONJA, CONJB>
    <<<dimGrid, dimBlock, 0, stream>>>(m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, offsetA, offsetB);
}
////////////////////////////////////////////////////////////////////////////////////////////////////
// NT, NC
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K, const int dim_vec,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int CONJA, const int CONJB>
void gemm_template_vbatched_nt(
    pastix_int_t* m, pastix_int_t* n, pastix_int_t* k,
    T const * const * dA_array, pastix_int_t* ldda,
    T const * const * dB_array, pastix_int_t* lddb,
    T**       dC_array, pastix_int_t* lddc,
    T alpha, T beta,
    pastix_int_t offsetA, pastix_int_t offsetB,
    pastix_int_t batchCount, cudaStream_t stream,
    pastix_int_t max_m, pastix_int_t max_n)
{
    dim3 dimBlock(DIM_X, DIM_Y);
    dim3 dimGrid( pastix_ceildiv( max_m, BLK_M ), pastix_ceildiv( max_n, BLK_N ), batchCount );
    gemm_template_vbatched_nt_kernel<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB, CONJA, CONJB>
    <<<dimGrid, dimBlock, 0, stream>>>(m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, offsetA, offsetB);
}
////////////////////////////////////////////////////////////////////////////////////////////////////
// TN, CN
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K, const int dim_vec,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int CONJA, const int CONJB>
void gemm_template_vbatched_tn(
    pastix_int_t* m, pastix_int_t* n, pastix_int_t* k,
    T const * const * dA_array, pastix_int_t* ldda,
    T const * const * dB_array, pastix_int_t* lddb,
    T**       dC_array, pastix_int_t* lddc,
    T alpha, T beta,
    pastix_int_t offsetA, pastix_int_t offsetB,
    pastix_int_t batchCount, cudaStream_t stream,
    pastix_int_t max_m, pastix_int_t max_n)
{
    dim3 dimBlock(DIM_X, DIM_Y);
    dim3 dimGrid( pastix_ceildiv( max_m, BLK_M ), pastix_ceildiv( max_n, BLK_N ), batchCount );
    gemm_template_vbatched_tn_kernel<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB, CONJA, CONJB>
    <<<dimGrid, dimBlock, 0, stream>>>(m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, offsetA, offsetB);
}
////////////////////////////////////////////////////////////////////////////////////////////////////
// TT, TC, CT, CC
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K, const int dim_vec,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int CONJA, const int CONJB>
void gemm_template_vbatched_tt(
    pastix_int_t* m, pastix_int_t* n, pastix_int_t* k,
    T const * const * dA_array, pastix_int_t* ldda,
    T const * const * dB_array, pastix_int_t* lddb,
    T**       dC_array, pastix_int_t* lddc,
    T alpha, T beta,
    pastix_int_t offsetA, pastix_int_t offsetB,
    pastix_int_t batchCount, cudaStream_t stream,
    pastix_int_t max_m, pastix_int_t max_n)
{
    dim3 dimBlock(DIM_X, DIM_Y);
    dim3 dimGrid( pastix_ceildiv( max_m, BLK_M ), pastix_ceildiv( max_n, BLK_N ), batchCount );
    gemm_template_vbatched_tt_kernel<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB, CONJA, CONJB>
    <<<dimGrid, dimBlock, 0, stream>>>(m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, offsetA, offsetB);
}


#endif //GEMM_TEMPLATE_KERNEL_VBATCHED_CUH
