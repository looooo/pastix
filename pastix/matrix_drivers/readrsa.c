/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "drivers.h"

/*
  Function: FORTRAN_CALL(wreadmtc)

  Declaration of the wreadmtc fortran function
  defined in skitf.f.

  Parameters:
     tmp1      - Maximum number of column                                    (INPUT)
     tmp2      - Maximum number of non zeros                                 (INPUT)
     tmp3      - job to be done (see skitf file)                             (INPUT)
     filename  - Path to the file to read from                               (INPUT)
     len       - length of *filname*                                         (INPUT)
     val       - Values of the elements of the matrix.                       (OUTPUT)
     row       - Rows of the elements of the matrix.                         (OUTPUT)
     col       - Index of first element of each column in *row* and *val*    (OUTPUT)
     crhs      - Right hand side(s).                                         (OUTPUT)
     nrhs      - Number of right hand side(s).                               (OUTPUT)
     RhsType   - Right hand side type                                        (OUTPUT)
     tmpNrow   - Number of rows.                                             (OUTPUT)
     tmpNcol   - Number of columns.                                          (OUTPUT)
     tmpNnzero - Number of non zeros.                                        (OUTPUT)
     title     - name of the matrix.                                         (OUTPUT)
     key       - key of the matrix (see skitf.f)                             (OUTPUT)
     Type      - Type of the matrix                                          (OUTPUT)
     ierr      - Error return value                                          (OUTPUT)

 */
void
FC_GLOBAL(wreadmtc,WREADMTC)(int        *tmp1,
                             int        *tmp2,
                             int        *tmp3,
                             const char *filename,
                             int        *len,
                             double     *val,
                             int        *row,
                             int        *col,
                             double     *crhs,
                             int        *nrhs,
                             char       *RhsType,
                             int        *tmpNrow,
                             int        *tmpNcol,
                             int        *tmpNnzero,
                             char       *title,
                             char       *key,
                             char       *Type,
                             int        *ierr);

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * readRSAHeader - Read the header structure of a RSA file
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The filename that contains the matrix stored in RSA format. This is
 *          an interface to the wreadmtx fortran fucntion provided within
 *          SparseKit.
 *
 * @param[out] N
 *          At exit, contains the number of rows/columns of the matrix.
 *
 * @param[out] Nnzero
 *          At exit, contains the number of non zero entries of the matrix.
 *
 * @param[out] Type
 *          At exit, contains the type of the matrix.
 *
 * @param[out] RhsType
 *          At exit, contains the type of the right hand side.
 *
 *******************************************************************************/
void readRSAHeader( const char *filename,
                    int        *N,
                    int        *Nnz,
                    char       *Type,
                    char       *RhsType )
{
    int     tmp;
    int    *col = NULL;
    int    *row = NULL;
    char    title[72+1];
    char    key[8+1];
    int     nrhs;
    int     len;
    int     ierr;
    double *val  = NULL;
    double *crhs = NULL;
    int     M;

    len = strlen(filename);
    tmp = 0;

    FC_GLOBAL(wreadmtc,WREADMTC)
        (&tmp, &tmp, &tmp, filename, &len, val, row, col, crhs, &nrhs,
         RhsType, &M, N, Nnz, title, key, Type, &ierr );

    if(ierr != 0) {
        fprintf(stderr, "cannot read matrix (job=0)\n");
    }

    if ( M != (*N) )
    {
        fprintf(stderr,"ERROR : (M != N)\n");
        exit(EXIT_FAILURE);
    }

    Type[3] = '\0';
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * readRSA - Read a RSA matrix file.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The filename that contains the matrix stored in RSA format. This is
 *          an interface to the wreadmtx fortran fucntion provided within
 *          SparseKit.
 *
 * @param[out] Nrow
 *          At exit, contains the number of rows of the matrix.
 *
 * @param[out] Ncol
 *          At exit, contains the number of columns of the matrix.
 *
 * @param[out] Nnzero
 *          At exit, contains the number of non zero entries of the matrix.
 *
 * @param[out] Type
 *          At exit, contains the type of the matrix.
 *
 * @param[out] RhsType
 *          At exit, contains the type of the right hand side.
 *
 *******************************************************************************/
void readRSA( const char   *filename,
              pastix_csc_t *csc )
{
    char    Type[4];
    char    RhsType[4];
    int     i;
    int     tmp;
    char    title[72+1];
    char    key[8+1];
    int     nrhs;
    int     len;
    int     ierr;
    double *crhs=NULL;
    int     M, N, Nnz;
    int    *tmpcolptr;
    int    *tmprows;
    int     base;

    readRSAHeader(filename, &N, &Nnz, Type, RhsType );

    tmpcolptr = (int*) malloc( (N+1) * sizeof(int) );
    assert( tmpcolptr );

    tmprows = (int*) malloc( Nnz * sizeof(int) );
    assert( tmprows );

    /* RSA reads only double */
    csc->flttype = PastixDouble;
    csc->avals = (double*) malloc( Nnz * sizeof(double) );
    assert( csc->avals );

    len  = strlen(filename);
    tmp  = 2;
    nrhs = 0;

    FC_GLOBAL(wreadmtc,WREADMTC)
        (&N, &Nnz, &tmp, filename, &len, csc->avals, tmprows, tmpcolptr, crhs,
         &nrhs, RhsType, &M, &N, &Nnz, title, key, Type, &ierr );

    base = (tmpcolptr[0] == 0) ? 0 : 1;
    assert( (tmpcolptr[N]-base) == Nnz );

    csc->gN = N;
    csc->n  = N;

    csc->colptr = (pastix_int_t*)malloc( (N+1) * sizeof(pastix_int_t) );
    assert( csc->colptr );
    for (i=0; i<N+1; i++)
        (csc->colptr)[i] = (pastix_int_t)(tmpcolptr[i]);
    memFree_null(tmpcolptr);

    csc->rows = (pastix_int_t*)malloc( Nnz * sizeof(pastix_int_t) );
    assert( csc->rows );
    for (i=0; i<Nnz; i++)
        (csc->rows)[i] = (pastix_int_t)(tmprows[i]);
    memFree_null(tmprows);

    RhsType[0]='\0';
    if(ierr != 0) {
        fprintf(stderr, "cannot read matrix (job=2)\n");
    }

    switch( Type[1] ){
    case 'S':
    case 's':
        csc->mtxtype = PastixSymmetric;
        break;
    case 'H':
    case 'h':
        csc->mtxtype = PastixHermitian;
        break;
    case 'U':
    case 'u':
    default:
        csc->mtxtype = PastixGeneral;
    }
    csc->flttype = PastixDouble;
}
