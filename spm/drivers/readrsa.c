/**
 * @file readrsa.c
 *
 *  $COPYRIGHTS$
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#include "common.h"
#include "spm_drivers.h"

/**
 * Function: FORTRAN_CALL(wreadmtc)
 *
 * Declaration of the wreadmtc fortran function
 * defined in skitf.f.
 *
 * Parameters:
 *    tmp1      - Maximum number of column                                    (INPUT)
 *    tmp2      - Maximum number of non zeros                                 (INPUT)
 *    tmp3      - job to be done (see skitf file)                             (INPUT)
 *    filename  - Path to the file to read from                               (INPUT)
 *    len       - length of *filname*                                         (INPUT)
 *    val       - Values of the elements of the matrix.                       (OUTPUT)
 *    row       - Rows of the elements of the matrix.                         (OUTPUT)
 *    col       - Index of first element of each column in *row* and *val*    (OUTPUT)
 *    crhs      - Right hand side(s).                                         (OUTPUT)
 *    nrhs      - Number of right hand side(s).                               (OUTPUT)
 *    RhsType   - Right hand side type                                        (OUTPUT)
 *    tmpNrow   - Number of rows.                                             (OUTPUT)
 *    tmpNcol   - Number of columns.                                          (OUTPUT)
 *    tmpNnzero - Number of non zeros.                                        (OUTPUT)
 *    title     - name of the matrix.                                         (OUTPUT)
 *    key       - key of the matrix (see skitf.f)                             (OUTPUT)
 *    Type      - Type of the matrix                                          (OUTPUT)
 *    ierr      - Error return value                                          (OUTPUT)
 *
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
 * @ingroup pastix_spm_driver
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
void
readRSAHeader( const char *filename,
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
 * @ingroup pastix_spm_driver
 *
 * readRSA - Read a RSA matrix file. This driver reads only real matrices, and
 * does not support complex matrices.
 * The matrix is returned in double, convert it to real if needed through TODO
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The filename that contains the matrix stored in RSA format. This is
 *          an interface to the wreadmtx fortran function provided within
 *          SparseKit.
 *
 * @param[in] spm
 *          At exit, contains the matrix in spm format.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the matrix has been read successfully
 *      \retval PASTIX_ERR_IO if a problem occured in the RSA driver
 *      \retval PASTIX_ERR_BADPARAMETER if the matrix is no in a supported format
 *
 *******************************************************************************/
int
readRSA( const char   *filename,
         pastix_spm_t *spm )
{
    char    Type[4];
    char    RhsType[4];
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

    readRSAHeader(filename, &N, &Nnz, Type, RhsType );

    switch( Type[1] ){
    case 'S':
    case 's':
        spm->mtxtype = PastixSymmetric;
        break;
    case 'H':
    case 'h':
        spm->mtxtype = PastixHermitian;
        /**
         * We should not arrive here, since the fortran driver is not able to
         * read complex matrices
         */
        fprintf(stderr,"readrsa: Unsupported Complex.\n");
        return PASTIX_ERR_BADPARAMETER;
    case 'U':
    case 'u':
        spm->mtxtype = PastixGeneral;
        break;
    default:
        fprintf(stderr,"readrsa: Unsupported type of matrix.\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    spm->flttype = PastixDouble;
    spm->fmttype = PastixCSC;
    spm->gN      = N;
    spm->n       = N;
    spm->gnnz    = Nnz;
    spm->nnz     = Nnz;
    spm->dof     = 1;
    spm->loc2glob= NULL;

    tmpcolptr = (int*) malloc( (N+1) * sizeof(int) );
    assert( tmpcolptr );

    tmprows = (int*) malloc( Nnz * sizeof(int) );
    assert( tmprows );

    /* RSA reads only double */
    spm->values = (double*) malloc( Nnz * sizeof(double) );
    assert( spm->values );

    len  = strlen(filename);
    tmp  = 2;
    nrhs = 0;

    FC_GLOBAL(wreadmtc,WREADMTC)
        (&N, &Nnz, &tmp, filename, &len, spm->values, tmprows, tmpcolptr, crhs,
         &nrhs, RhsType, &M, &N, &Nnz, title, key, Type, &ierr );

    assert( (tmpcolptr[N]-tmpcolptr[0]) == Nnz );

    spm->colptr  = spmIntConvert( N+1, tmpcolptr );
    spm->rowptr  = spmIntConvert( Nnz, tmprows );

    RhsType[0] = '\0';
    if(ierr != 0) {
        fprintf(stderr, "cannot read matrix (job=2)\n");
        return PASTIX_ERR_IO;
    }

    return PASTIX_SUCCESS;
}
