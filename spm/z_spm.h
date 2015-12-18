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

int    z_spmGeCSCv(int trans, pastix_complex64_t alpha, const pastix_csc_t *csc, const pastix_complex64_t *x, pastix_complex64_t beta, pastix_complex64_t *y);
int    z_spmSyCSCv(           pastix_complex64_t alpha, const pastix_csc_t *csc, const pastix_complex64_t *x, pastix_complex64_t beta, pastix_complex64_t *y);
int    z_spmHeCSCv(           pastix_complex64_t alpha, const pastix_csc_t *csc, const pastix_complex64_t *x, pastix_complex64_t beta, pastix_complex64_t *y);

double       z_spmNorm( int ntype, const pastix_csc_t *csc );
void         z_spmSort( pastix_csc_t *csc );
pastix_int_t z_spmMergeDuplicate( pastix_csc_t *csc );
pastix_int_t z_spmSymmetrize( pastix_csc_t *csc );

int z_spmGenRHS(int type, int nrhs, const pastix_csc_t *spm, void *x, int ldx, void *b, int ldb );
int z_spmCheckAxb( int nrhs, const pastix_csc_t *spm, void *x0, int ldx0, void *b, int ldb, const void *x, int ldx );

#endif /* _z_spm_H_ */
