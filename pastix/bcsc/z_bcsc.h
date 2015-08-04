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

int z_bcscHemv(      pastix_complex64_t  alpha,
               const pastix_bcsc_t      *bcsc,
               const pastix_complex64_t *x,
                     pastix_complex64_t  beta,
                     pastix_complex64_t *y );

double z_bcscNorm( pastix_normtype_t ntype, const pastix_bcsc_t *bcsc );

double z_bcscBerr( void         *r1,
                   void         *r2,
                   pastix_int_t  n );

double z_bcscNormErr( void         *r1,
                      void         *r2,
                      pastix_int_t  n );

int z_bcscScal( void               *x,
                pastix_complex64_t  alpha,
                pastix_int_t        n,
                pastix_int_t        smxnbr );

int z_bcscAxpy(pastix_complex64_t  alpha,
               void               *x,
               pastix_int_t        n,
               void               *y,
               pastix_int_t        smxnbr );

double z_bcscAxpb( const pastix_bcsc_t *bcsc,
                   void                *x,
                   void                *b );

double z_bcscDotc( void                *x,
                   void                *y,
                   pastix_int_t         n);

double z_vectFrobeniusNorm( void *, pastix_int_t );

int z_bcscApplyPerm( pastix_int_t m,
                     pastix_int_t n,
                     pastix_complex64_t *A,
                     pastix_int_t lda,
                     pastix_int_t *perm );
#endif /* _Z_BCSC_H_ */
