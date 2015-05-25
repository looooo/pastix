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

int spmConvert( int ofmttype, pastix_csc_t *ospm );

int genRHS(pastix_csc_t *csc, void **rhs );

pastix_int_t spmFindBase( const pastix_csc_t *spm );
double       spmNorm( int ntype, const pastix_csc_t *csc );
void         spmExit( pastix_csc_t *spm );

#endif /* _CSC_H_ */
