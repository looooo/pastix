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
    pastix_int_t *rows;      /*< List of edges for each vertex                       */
    pastix_int_t *loc2glob;  /*< Corresponding numbering from local to global        */
    void         *avals;     /*< Values stored in the matrix                         */
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

/* Not implemented */
int z_spmConvertCSC2CSR( int ofmttype, pastix_csc_t *spm );
int c_spmConvertCSC2CSR( int ofmttype, pastix_csc_t *spm );
int d_spmConvertCSC2CSR( int ofmttype, pastix_csc_t *spm );
int s_spmConvertCSC2CSR( int ofmttype, pastix_csc_t *spm );
int p_spmConvertCSC2CSR( int ofmttype, pastix_csc_t *spm );

/* Not implemented */
int z_spmConvertCSC2IJV( int ofmttype, pastix_csc_t *spm );
int c_spmConvertCSC2IJV( int ofmttype, pastix_csc_t *spm );
int d_spmConvertCSC2IJV( int ofmttype, pastix_csc_t *spm );
int s_spmConvertCSC2IJV( int ofmttype, pastix_csc_t *spm );
int p_spmConvertCSC2IJV( int ofmttype, pastix_csc_t *spm );

/* Not implemented */
int z_spmConvertCSR2CSC( int ofmttype, pastix_csc_t *spm );
int c_spmConvertCSR2CSC( int ofmttype, pastix_csc_t *spm );
int d_spmConvertCSR2CSC( int ofmttype, pastix_csc_t *spm );
int s_spmConvertCSR2CSC( int ofmttype, pastix_csc_t *spm );
int p_spmConvertCSR2CSC( int ofmttype, pastix_csc_t *spm );

/* Not implemented */
int z_spmConvertCSR2IJV( int ofmttype, pastix_csc_t *spm );
int c_spmConvertCSR2IJV( int ofmttype, pastix_csc_t *spm );
int d_spmConvertCSR2IJV( int ofmttype, pastix_csc_t *spm );
int s_spmConvertCSR2IJV( int ofmttype, pastix_csc_t *spm );
int p_spmConvertCSR2IJV( int ofmttype, pastix_csc_t *spm );

int z_spmConvertIJV2CSC( int ofmttype, pastix_csc_t *spm );
int c_spmConvertIJV2CSC( int ofmttype, pastix_csc_t *spm );
int d_spmConvertIJV2CSC( int ofmttype, pastix_csc_t *spm );
int s_spmConvertIJV2CSC( int ofmttype, pastix_csc_t *spm );
int p_spmConvertIJV2CSC( int ofmttype, pastix_csc_t *spm );

/* Not implemented */
int z_spmConvertIJV2CSR( int ofmttype, pastix_csc_t *spm );
int c_spmConvertIJV2CSR( int ofmttype, pastix_csc_t *spm );
int d_spmConvertIJV2CSR( int ofmttype, pastix_csc_t *spm );
int s_spmConvertIJV2CSR( int ofmttype, pastix_csc_t *spm );
int p_spmConvertIJV2CSR( int ofmttype, pastix_csc_t *spm );


int spmConvert( int ofmttype, pastix_csc_t *ospm );

#endif /* _CSC_H_ */
