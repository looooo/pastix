/**
 * @file pastix_zcuda.h
 *
 * PaStiX GPU kernel header.
 *
 * @copyright 2011-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2021-08-25
 * @precisions normal z -> c d s
 *
 */
#ifndef _pastix_zcuda_h_
#define _pastix_zcuda_h_

#if defined(PASTIX_WITH_CUDA)

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuComplex.h>

/**
 * @addtogroup kernel_blas_lapack
 * @{
 *    @name PastixComplex64 cblk-BLAS GPU kernels
 *    @{
 */
pastix_fixdbl_t gpucblk_zgemmsp( pastix_coefside_t sideA, pastix_trans_t trans,
                                 const SolverCblk *cblk, const SolverBlok *blok, SolverCblk *fcblk,
                                 const cuDoubleComplex *A, const cuDoubleComplex *B, cuDoubleComplex *C,
                                 const pastix_lr_t *lowrank, cudaStream_t stream );

pastix_fixdbl_t gpublok_zgemmsp( pastix_trans_t trans,
                                 const SolverCblk *cblk, SolverCblk *fcblk,
                                 pastix_int_t blok_mk, pastix_int_t blok_nk, pastix_int_t blok_mn,
                                 const void *A, const void *B, void *C,
                                 const pastix_lr_t *lowrank, cudaStream_t stream );

pastix_fixdbl_t gpublok_ztrsmsp( pastix_side_t side, pastix_uplo_t uplo,
                                 pastix_trans_t trans, pastix_diag_t diag,
                                 const SolverCblk *cblk, pastix_int_t blok_m,
                                 const void *A, void *C,
                                 const pastix_lr_t *lowrank, cudaStream_t stream );

pastix_fixdbl_t gpu_zgemmsp_fermi( const SolverMatrix *solvmatr,
                                   pastix_uplo_t uplo, pastix_trans_t trans,
                                   int *blocktab,
                                   const SolverCblk      *cblk,
                                   const SolverBlok      *blok,
                                   SolverCblk      *fcblk,
                                   const cuDoubleComplex *A,
                                   const cuDoubleComplex *B,
                                   cuDoubleComplex *C,
                                   cudaStream_t stream );

/**
 *    @}
 * @}
 */
#endif

#endif /* _pastix_zcuda_h_ */
