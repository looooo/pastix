/**
 * @file readhb.c
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
#include <stdio.h>
#include "common.h"
#include "spm_drivers.h"
#include "drivers/iohb.h"

/**
 * ******************************************************************************
 *
 * @ingroup pastix_spm_driver
 *
 * readHB - Interface to the Harwell-Boeing C driver (iohb.c)
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The file containing the matrix.
 *
 * @param[in] spm
 *          At exit, contains the matrix in spm format.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the matrix has been read successfully
 *      \retval PASTIX_ERR_IO if a problem occured in the Harwell Boeing driver
 *      \retval PASTIX_ERR_BADPARAMETER if the matrix is no in a supported format
 *
 *******************************************************************************/
int
readHB( const char   *filename,
        pastix_spm_t *spm )
{
    int M, N, nz, nrhs;

    /* Harwell Boeing is a variant of RSA */
    spm->fmttype = PastixCSC;
    spm->dof     = 1;
    spm->loc2glob= NULL;

    /* Read header informations */
    {
        char *Type;
        Type = malloc(4*sizeof(char));
        Type[0] ='a';

        readHB_info(filename, &M, &N, &nz, &Type, &nrhs);

        if ( M != N ) {
            fprintf(stderr, "readHB: PaStiX does not support non square matrices (m=%d, N=%d\n", M, N);
            return PASTIX_ERR_BADPARAMETER;
        }

        spm->gN   = M;
        spm->n    = M;
        spm->gnnz = nz;
        spm->nnz  = nz;

        /* Check float type */
        switch( Type[0] ) {
        case 'C':
        case 'c':
            spm->flttype = PastixComplex64;
            break;
        case 'R':
        case 'r':
            spm->flttype = PastixDouble;
            break;
        case 'P':
        case 'p':
            spm->flttype = PastixPattern;
            break;
        default:
            fprintf(stderr, "readhb: Floating type unknown (%c)\n", Type[0]);
            return PASTIX_ERR_BADPARAMETER;
        }

        /* Check Symmetry */
        switch( Type[1] ) {
        case 'S':
        case 's':
            spm->mtxtype = PastixSymmetric;
            break;
        case 'H':
        case 'h':
            spm->mtxtype = PastixHermitian;
            assert( spm->flttype == PastixDouble );
            break;
        case 'U':
        case 'u':
        default:
            spm->mtxtype = PastixGeneral;
        }
        free(Type);
    }

    /* Read the matrix and its values */
    {
        int    *colptr, *rowind;
        int     rc;

        rc = readHB_newmat_double( filename, &M, &N, &nz,
                                   &colptr, &rowind, (double**)(&(spm->values)) );

        if (rc == 0) {
            fprintf(stderr, "readhb: Error in reading the HB matrix values\n");
            return PASTIX_ERR_IO;
        }

        /* Move the colptr/rowind from int to pastix_int_t if different sizes */
        spm->colptr = spmIntConvert(spm->n+1, colptr);
        spm->rowptr = spmIntConvert(spm->nnz, rowind);
    }
    return PASTIX_SUCCESS;
}
