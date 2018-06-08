/**
 * @file bcsc_z.h
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 * @precisions normal z -> c d s
 *
 **/
#ifndef _bcsc_z_h_
#define _bcsc_z_h_

/**
 * @addtogroup bcsc_internal
 * @{
 *
 *    @name PastixComplex64 initialization functions
 *    @{
 */
void bcsc_zinit_centralized( const spmatrix_t     *spm,
                             const pastix_order_t *ord,
                             const SolverMatrix   *solvmtx,
                             const pastix_int_t   *col2cblk,
                                   int             initAt,
                                   pastix_bcsc_t  *bcsc );

void bcsc_zsort( const pastix_bcsc_t *bcsc,
                 pastix_int_t        *rowtab,
                 pastix_complex64_t  *valtab );

/**
 *   @}
 * @}
 *
 * @addtogroup bcsc
 * @{
 *
 *    @name PastixComplex64 vector(s) operations
 *    @{
 */
double bvec_znrm2( pastix_int_t              n,
                   const pastix_complex64_t *x );
int    bvec_zscal( pastix_int_t        n,
                   pastix_complex64_t  alpha,
                   pastix_complex64_t *x );
int    bvec_zaxpy( pastix_int_t              n,
                   pastix_complex64_t        alpha,
                   const pastix_complex64_t *x,
                   pastix_complex64_t       *y );
#if defined(PRECISION_z) || defined(PRECISION_c)
pastix_complex64_t bvec_zdotc( pastix_int_t              n,
                               const pastix_complex64_t *x,
                               const pastix_complex64_t *y );
#endif
pastix_complex64_t bvec_zdotu( pastix_int_t              n,
                               const pastix_complex64_t *x,
                               const pastix_complex64_t *y );

int bvec_zswap( pastix_int_t m,
                pastix_int_t n,
                pastix_complex64_t *A,
                pastix_int_t lda,
                pastix_int_t *perm );

/**
 *    @}
 *
 *    @name PastixComplex64 matrix operations
 *    @{
 */
double bcsc_znorm( pastix_normtype_t    ntype,
                   const pastix_bcsc_t *bcsc );

int    bcsc_zspmv(       pastix_trans_t      trans,
                         pastix_complex64_t  alpha,
                   const pastix_bcsc_t      *bcsc,
                   const pastix_complex64_t *x,
                         pastix_complex64_t  beta,
                         pastix_complex64_t *y );
/**
 *    @}
 * @}
 */
#endif /* _bcsc_z_h_ */
