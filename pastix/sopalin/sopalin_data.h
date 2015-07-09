/**
 *
 * @file sopalin_data.h
 *
 *  PaStiX sopalin routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _SOPALIN_DATA_H_
#define _SOPALIN_DATA_H_

struct sopalin_data_s {
    void         *sched;
    SolverMatrix *solvmtx;
    double        diagthreshold; /* Threshold for static pivoting on diagonal value */
};
typedef struct sopalin_data_s sopalin_data_t;

void coeftab_zdump( const SolverMatrix *solvmtx, const char *filename );
void coeftab_cdump( const SolverMatrix *solvmtx, const char *filename );
void coeftab_ddump( const SolverMatrix *solvmtx, const char *filename );
void coeftab_sdump( const SolverMatrix *solvmtx, const char *filename );

void sequential_ztrsm( int side, int uplo, int trans, int diag, sopalin_data_t *sopalin_data, int nrhs, pastix_complex64_t *b, int ldb );
void sequential_ctrsm( int side, int uplo, int trans, int diag, sopalin_data_t *sopalin_data, int nrhs, pastix_complex32_t *b, int ldb );
void sequential_dtrsm( int side, int uplo, int trans, int diag, sopalin_data_t *sopalin_data, int nrhs, double *b, int ldb );
void sequential_strsm( int side, int uplo, int trans, int diag, sopalin_data_t *sopalin_data, int nrhs, float *b, int ldb );

void sequential_zdiag( sopalin_data_t *sopalin_data, int nrhs, pastix_complex64_t *b, int ldb );
void sequential_cdiag( sopalin_data_t *sopalin_data, int nrhs, pastix_complex32_t *b, int ldb );
void sequential_ddiag( sopalin_data_t *sopalin_data, int nrhs, double *b, int ldb );
void sequential_sdiag( sopalin_data_t *sopalin_data, int nrhs, float *b, int ldb );

void sopalin_zgetrf( sopalin_data_t *sopalin_data );
void sopalin_cgetrf( sopalin_data_t *sopalin_data );
void sopalin_dgetrf( sopalin_data_t *sopalin_data );
void sopalin_sgetrf( sopalin_data_t *sopalin_data );

void sopalin_zhetrf( sopalin_data_t *sopalin_data );
void sopalin_chetrf( sopalin_data_t *sopalin_data );

void sopalin_zpotrf( sopalin_data_t *sopalin_data );
void sopalin_cpotrf( sopalin_data_t *sopalin_data );
void sopalin_dpotrf( sopalin_data_t *sopalin_data );
void sopalin_spotrf( sopalin_data_t *sopalin_data );

void sopalin_zsytrf( sopalin_data_t *sopalin_data );
void sopalin_csytrf( sopalin_data_t *sopalin_data );
void sopalin_dsytrf( sopalin_data_t *sopalin_data );
void sopalin_ssytrf( sopalin_data_t *sopalin_data );

#endif /* _SOPALIN_DATA_H_ */


