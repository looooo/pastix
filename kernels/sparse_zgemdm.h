/**
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#ifndef SPARSE_ZGEMDM_GPU_H
#define SPARSE_ZGEMDM_GPU_H

#define sparse_zgemdm_blocktab_kernel_N_T_64_16_4_16_4 \
  PASTIX_PREFIX_F(sparse_zgemdm_blocktab_kernel_N_T_64_16_4_16_4)
#define magmablas_sparse_zgemdm_blocktab_kernel_N_T_64_16_4_16_4 \
  PASTIX_PREFIX_F(magmablas_sparse_zgemdm_blocktab_kernel_N_T_64_16_4_16_4)


#ifdef __cplusplus
extern "C"
__global__
#endif
/*
 * Function: sparse_zgemdm_kernel_N_T_64_16_4_16_4
 *
 * CUDA kernel to update C.
 *
 * Performs : $C \leftarrow \alpha A \times D \times B^T + \beta B$
 *
 * Parameters:
 *   m          - Number of rows in *C*.
 *   n          - Number of columns in *C*.
 *   k          - Number of columns in *A*.
 *   alpha      - A coefficient .
 *   A          - $m \times k$ matrix.
 *   lda        - Leading dimension of *A*.
 *   D          - Diagonal of the column from which *A* is a part.
 *   ldd        - Stride between two entries of *D*.
 *   B          - $n \times k$ matrix.
 *   ldb        - Leading dimension of *B*.
 *   beta       - A coefficient.
 *   C          - $m \times n$ matrix.
 *   ldc        - Leading dimension of *C*.
 *   blocknbr   - Number of blocks in *A*.
 *   blocktab   - Array containing first and last row of each block in *A*.
 *                blocktab[2i]   : first row of block i,
 *                blocktab[2i+1] : last row of block i.
 *   fblocknbr  - number of blocks in *C*.
 *   fblocktab  - Array containing first and last row of each block in *C*.
 */
void
sparse_zgemdm_kernel_N_T_64_16_4_16_4(int m, int n, int k,
                                     cuDoubleComplex alpha,
                                     const cuDoubleComplex *A, int lda,
                                     const cuDoubleComplex *D, int ldd,
                                     const cuDoubleComplex *B, int ldb,
                                     cuDoubleComplex beta,
                                     cuDoubleComplex       *C, int ldc,
                                     int blocknbr,  const int * blocktab,
                                     int fblocknbr, const int * fblocktab);
#ifdef __cplusplus
extern "C"
#endif

/*
 * Function: magmablas_sparse_zgemdm_kernel_N_T_64_16_4_16_4
 *
 * Interface to the CUDA kernel <sparse_zgemdm_kernel_N_T_64_16_4_16_4>.
 *
 * Parameters:
 *   m          - Number of rows in *C*.
 *   n          - Number of columns in *C*.
 *   k          - Number of columns in *A*.
 *   alpha      - A coefficient .
 *   A          - $m \times k$ matrix.
 *   lda        - Leading dimension of *A*.
 *   B          - $n \times k$ matrix.
 *   ldb        - Leading dimension of *B*.
 *   D          - Diagonal of the column from which *A* is a part.
 *   ldd        - Stride between two entries of *D*.
 *   beta       - A coefficient.
 *   C          - $m \times n$ matrix.
 *   ldc        - Leading dimension of *C*.
 *   blocknbr   - Number of blocks in *A*.
 *   blocktab   - Array containing first and last row of each block in *A*.
 *                blocktab[2i]   : first row of block i,
 *                blocktab[2i+1] : last row of block i.
 *   fblocknbr  - number of blocks in *C*.
 *   fblocktab  - Array containing first and last row of each block in *C*.
 */
void
magmablas_sparse_zgemdm_kernel_N_T_64_16_4_16_4(int m, int n, int k,
                                                cuDoubleComplex alpha,
                                                const cuDoubleComplex *A, int lda,
                                                const cuDoubleComplex *D, int ldd,
                                                const cuDoubleComplex *B, int ldb,
                                                cuDoubleComplex beta,
                                                cuDoubleComplex       *c, int ldc,
                                                int blocknbr,  const int * blocktab,
                                                int fblocknbr, const int * fblocktab);
#endif
