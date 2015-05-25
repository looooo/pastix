/**
 *
 * @file z_spm.h
 *
 *  PaStiX sparse matrix routines to handle different format of sparse matrices.
 *  $COPYRIGHTS$
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> c d s p
 *
 **/
#ifndef _z_spm_H_
#define _z_spm_H_

int z_spmConvertCSC2CSR( pastix_csc_t *spm );
int z_spmConvertCSC2IJV( pastix_csc_t *spm );
int z_spmConvertCSR2CSC( pastix_csc_t *spm );
int z_spmConvertCSR2IJV( pastix_csc_t *spm );
int z_spmConvertIJV2CSC( pastix_csc_t *spm );
int z_spmConvertIJV2CSR( pastix_csc_t *spm );

int    z_spmGeCSCv(char trans, pastix_complex64_t alpha, pastix_csc_t *csc, pastix_complex64_t *x, pastix_complex64_t beta, pastix_complex64_t *b);
int    z_spmSyCSCv(pastix_complex64_t alpha, pastix_csc_t *csc, pastix_complex64_t *x, pastix_complex64_t beta, pastix_complex64_t *b);
int    z_spmHeCSCv(pastix_complex64_t alpha, pastix_csc_t *csc, pastix_complex64_t *x, pastix_complex64_t beta, pastix_complex64_t *b);
int    z_spm_genRHS(pastix_csc_t *csc, void **rhs );
double z_spmNorm( int ntype, const pastix_csc_t *csc );

#endif /* _z_spm_H_ */
