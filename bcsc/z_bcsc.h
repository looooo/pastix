/**
 * @file z_bcsc.h
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
#ifndef _z_bcsc_h_
#define _z_bcsc_h_

void bcsc_zinit_centralized( const spmatrix_t     *spm,
                             const pastix_order_t *ord,
                             const SolverMatrix   *solvmtx,
                             const pastix_int_t   *col2cblk,
                                   int             initAt,
                                   pastix_bcsc_t  *bcsc );


double bcsc_znrm2( pastix_int_t              n,
                   const pastix_complex64_t *x );
int    bcsc_zscal( pastix_int_t        n,
                   pastix_complex64_t  alpha,
                   pastix_complex64_t *x );
int    bcsc_zaxpy( pastix_int_t              n,
                   pastix_complex64_t        alpha,
                   const pastix_complex64_t *x,
                   pastix_complex64_t       *y );
#if defined(PRECISION_z) || defined(PRECISION_c)
pastix_complex64_t bcsc_zdotc( pastix_int_t              n,
                               const pastix_complex64_t *x,
                               const pastix_complex64_t *y );
#endif
pastix_complex64_t bcsc_zdotu( pastix_int_t              n,
                               const pastix_complex64_t *x,
                               const pastix_complex64_t *y );




void z_bcscSort( const pastix_bcsc_t *bcsc,
                 pastix_int_t        *rowtab,
                 pastix_complex64_t  *valtab );

double z_bcscNorm( pastix_normtype_t ntype,
                   const pastix_bcsc_t *bcsc );

int z_bcscGemv(      pastix_trans_t      trans,
                     pastix_complex64_t  alpha,
               const pastix_bcsc_t      *bcsc,
               const pastix_complex64_t *x,
                     pastix_complex64_t  beta,
                     pastix_complex64_t *y );

double z_bcscBerr( pastix_complex64_t *r1,
                   pastix_complex64_t *r2,
                   pastix_int_t        n );

void z_bcscAxpb( pastix_trans_t       trans,
                 const pastix_bcsc_t *bcsc,
                 pastix_complex64_t  *x,
                 pastix_complex64_t  *b,
                 pastix_complex64_t  *r );

int z_bcscApplyPerm( pastix_int_t m,
                     pastix_int_t n,
                     pastix_complex64_t *A,
                     pastix_int_t lda,
                     pastix_int_t *perm );

#endif /* _z_bcsc_h_ */
