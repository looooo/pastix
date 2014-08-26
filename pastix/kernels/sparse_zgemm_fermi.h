/**
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/

#ifndef SPARSE_GEMM_FERMI_H
#define SPARSE_GEMM_FERMI_H

#if (!defined PRECISION_z && !defined PRECISION_c && !defined PRECISION_d && !defined PRECISION_s)
/* Required because CUDA compiler does not get the -DPRECISION_x */
#  define PRECISION_z
#endif
#include <cuda.h>
#if defined(PRECISION_z) || defined(PRECISION_c)
#include <cuComplex.h>
#endif
#include <cuda_runtime.h>

#define magmablas_zgemm_fermi  magmablas_zgemm
#define magmablas_zgemdm_fermi magmablas_zgemdm

////////////////////////////////////////////////////////////////////////////////
#ifdef PRECISION_z
#  define GENERATE_SM_VERSION_NAME_I(func, version) sparse_z##func##_SM##version
#endif
#ifdef PRECISION_c
#  define GENERATE_SM_VERSION_NAME_I(func, version) sparse_c##func##_SM##version
#  define cublasZgemm cublasCgemm
#  define cuDoubleComplex cuFloatComplex
#endif
#ifdef PRECISION_d
#  define GENERATE_SM_VERSION_NAME_I(func, version) sparse_d##func##_SM##version
#  define cublasZgemm cublasDgemm
#  define cuDoubleComplex double
#endif
#ifdef PRECISION_s
#  define GENERATE_SM_VERSION_NAME_I(func, version) sparse_s##func##_SM##version
#  define cublasZgemm cublasSgemm
#  define cuDoubleComplex float
#endif
#define GENERATE_SM_VERSION_KERNEL_NAME_I2(func, version) GENERATE_SM_VERSION_NAME_I(func, version)
#define GENERATE_SM_VERSION_KERNEL_NAME(func) GENERATE_SM_VERSION_NAME_I2(func, CUDA_SM_VERSION)


#define GENERATE_SM_VERSION_NAME_I2(func, version) GENERATE_SM_VERSION_NAME_I(func, version)
#define GENERATE_SM_VERSION_NAME(func) GENERATE_SM_VERSION_NAME_I2(func, CUDA_SM_VERSION)

#ifdef __cplusplus
extern "C"
#endif
void
GENERATE_SM_VERSION_NAME(gemm)( char TRANSA, char TRANSB,
                                int m , int n , int k ,
                                cuDoubleComplex alpha,
                                const cuDoubleComplex *d_A, int lda,
                                const cuDoubleComplex *d_B, int ldb,
                                cuDoubleComplex beta,
                                cuDoubleComplex *d_C, int ldc,
                                int blocknbr, const int *blocktab,
                                int fblocknbr, const int *fblocktab,
                                cudaStream_t stream );
#ifdef __cplusplus
extern "C"
#endif
void
GENERATE_SM_VERSION_NAME(gemdm)( char TRANSA, char TRANSB,
                                 int m , int n , int k ,
                                 cuDoubleComplex alpha,
                                 const cuDoubleComplex *d_A, int lda,
                                 const cuDoubleComplex *d_D, int ldd,
                                 const cuDoubleComplex *d_B, int ldb,
                                 cuDoubleComplex beta,
                                 cuDoubleComplex *d_C, int ldc,
                                 int blocknbr, const int *blocktab,
                                 int fblocknbr, const int *fblocktab,
                                 cudaStream_t stream );

#endif
