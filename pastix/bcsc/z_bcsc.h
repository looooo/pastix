/**
 * @file z_bcsc.h
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 * @precisions normal z -> c d s
 *
 **/
#ifndef _Z_BCSC_H_
#define _Z_BCSC_H_

void
z_bcscSort( const pastix_bcsc_t *bcsc,
            pastix_int_t        *rowtab,
            pastix_complex64_t  *valtab );

void z_bcscInitCentralized( const pastix_csc_t  *csc,
                            const Order         *ord,
                            const SolverMatrix  *solvmtx,
                            const pastix_int_t  *col2cblk,
                                  int            initAt,
                                  pastix_bcsc_t *bcsc );

int z_bcscGemv(      pastix_trans_t      trans,
                     pastix_complex64_t  alpha,
               const pastix_bcsc_t      *bcsc,
               const pastix_complex64_t *x,
                     pastix_complex64_t  beta,
                     pastix_complex64_t *y );

int z_bcscSymv(      pastix_complex64_t  alpha,
               const pastix_bcsc_t      *bcsc,
               const pastix_complex64_t *x,
                     pastix_complex64_t  beta,
                     pastix_complex64_t *y );

#if defined(PRECISION_z) || defined(PRECISION_c)
int z_bcscHemv(      pastix_complex64_t  alpha,
               const pastix_bcsc_t      *bcsc,
               const pastix_complex64_t *x,
                     pastix_complex64_t  beta,
                     pastix_complex64_t *y );
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

double z_bcscNorm( pastix_normtype_t ntype, const pastix_bcsc_t *bcsc );

int z_bcscApplyPerm( pastix_int_t m,
                     pastix_int_t n,
                     pastix_complex64_t *A,
                     pastix_int_t lda,
                     pastix_int_t *perm );
#endif /* _Z_BCSC_H_ */
