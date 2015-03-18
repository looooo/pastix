/**
 *
 * @file csc.h
 *
 *  PaStiX sparse matrix routines to handle different format of sparse matrices.
 *  $COPYRIGHTS$
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _CSC_H_
#define _CSC_H_

/**
 * @ingroup pastix_spm
 *
 * @struct pastix_spm_s - Sparse matrix data structure
 */
struct pastix_spm_s {
    int           mtxtype;   /*< Matrix structure: PastixGeneral, PastixSymmetric or PastixHermitian.          */
    int           flttype;   /*< avals datatype: PastixFloat, PastixDouble, PastixComplex32 or PastixComplex64 */
    int           fmttype;   /*< Matrix storage format: PastixCSC, PastixCSR, PastixIJV                        */
    pastix_int_t  gN;        /*< Global number of vertices in the compressed graph   */
    pastix_int_t  n;         /*< Local number of vertices in the compressed graph    */
    pastix_int_t  gnnz;      /*< Global number of non zeroes in the compressed graph */
    pastix_int_t  nnz;       /*< Local number of non zeroes in the compressed graph  */
    pastix_int_t  dof;       /*< Number of degrees of freedom per unknown            */
    pastix_int_t *colptr;    /*< List of indirections to rows for each vertex        */
    pastix_int_t *rowptr;    /*< List of edges for each vertex                       */
    pastix_int_t *loc2glob;  /*< Corresponding numbering from local to global        */
    void         *values;    /*< Values stored in the matrix                         */
};
typedef struct pastix_spm_s pastix_csc_t;

int
csc_load( pastix_int_t  *n,
          pastix_int_t **colptr,
          pastix_int_t **rows,
          int           *valtype,
          void         **values,
          int           *dof,
          FILE          *infile );

int
csc_save( pastix_int_t  n,
          pastix_int_t *colptr,
          pastix_int_t *rows,
          int           ft,
          void         *values,
          int           dof,
          FILE         *outfile );

int cscLoad( pastix_csc_t *csc, FILE *infile );
int cscSave( pastix_csc_t *csc, FILE *outfile );

int z_spmConvertCSC2CSR( pastix_csc_t *spm );
int c_spmConvertCSC2CSR( pastix_csc_t *spm );
int d_spmConvertCSC2CSR( pastix_csc_t *spm );
int s_spmConvertCSC2CSR( pastix_csc_t *spm );
int p_spmConvertCSC2CSR( pastix_csc_t *spm );

int z_spmConvertCSC2IJV( pastix_csc_t *spm );
int c_spmConvertCSC2IJV( pastix_csc_t *spm );
int d_spmConvertCSC2IJV( pastix_csc_t *spm );
int s_spmConvertCSC2IJV( pastix_csc_t *spm );
int p_spmConvertCSC2IJV( pastix_csc_t *spm );

int z_spmConvertCSR2CSC( pastix_csc_t *spm );
int c_spmConvertCSR2CSC( pastix_csc_t *spm );
int d_spmConvertCSR2CSC( pastix_csc_t *spm );
int s_spmConvertCSR2CSC( pastix_csc_t *spm );
int p_spmConvertCSR2CSC( pastix_csc_t *spm );

int z_spmConvertCSR2IJV( pastix_csc_t *spm );
int c_spmConvertCSR2IJV( pastix_csc_t *spm );
int d_spmConvertCSR2IJV( pastix_csc_t *spm );
int s_spmConvertCSR2IJV( pastix_csc_t *spm );
int p_spmConvertCSR2IJV( pastix_csc_t *spm );

int z_spmConvertIJV2CSC( pastix_csc_t *spm );
int c_spmConvertIJV2CSC( pastix_csc_t *spm );
int d_spmConvertIJV2CSC( pastix_csc_t *spm );
int s_spmConvertIJV2CSC( pastix_csc_t *spm );
int p_spmConvertIJV2CSC( pastix_csc_t *spm );

int z_spmConvertIJV2CSR( pastix_csc_t *spm );
int c_spmConvertIJV2CSR( pastix_csc_t *spm );
int d_spmConvertIJV2CSR( pastix_csc_t *spm );
int s_spmConvertIJV2CSR( pastix_csc_t *spm );
int p_spmConvertIJV2CSR( pastix_csc_t *spm );

int spmConvert( int ofmttype, pastix_csc_t *ospm );

int z_spmGeCSCv(char trans, pastix_complex64_t alpha, pastix_csc_t *csc, pastix_complex64_t *x, pastix_complex64_t beta, pastix_complex64_t *b);
int c_spmGeCSCv(char trans, pastix_complex32_t alpha, pastix_csc_t *csc, pastix_complex32_t *x, pastix_complex32_t beta, pastix_complex32_t *b);
int d_spmGeCSCv(char trans, double alpha, pastix_csc_t *csc, double *x, double beta, double *b);
int s_spmGeCSCv(char trans, float alpha, pastix_csc_t *csc, float *x, float beta, float *b);

int z_spmSyCSCv(pastix_complex64_t alpha, pastix_csc_t *csc, pastix_complex64_t *x, pastix_complex64_t beta, pastix_complex64_t *b);
int c_spmSyCSCv(pastix_complex32_t alpha, pastix_csc_t *csc, pastix_complex32_t *x, pastix_complex32_t beta, pastix_complex32_t *b);
int d_spmSyCSCv(double alpha, pastix_csc_t *csc, double *x, double beta, double *b);
int s_spmSyCSCv(float alpha, pastix_csc_t *csc, float *x, float beta, float *b);

int z_spmHeCSCv(pastix_complex64_t alpha, pastix_csc_t *csc, pastix_complex64_t *x, pastix_complex64_t beta, pastix_complex64_t *b);
int c_spmHeCSCv(pastix_complex32_t alpha, pastix_csc_t *csc, pastix_complex32_t *x, pastix_complex32_t beta, pastix_complex32_t *b);

int z_spm_genRHS(pastix_csc_t *csc, void **rhs );
int c_spm_genRHS(pastix_csc_t *csc, void **rhs );
int d_spm_genRHS(pastix_csc_t *csc, void **rhs );
int s_spm_genRHS(pastix_csc_t *csc, void **rhs );
int genRHS(pastix_csc_t *csc, void **rhs );

pastix_int_t spmFindBase( pastix_csc_t *spm );
#endif /* _CSC_H_ */
