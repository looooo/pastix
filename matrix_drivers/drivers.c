/**
 * @file drivers.c
 *
 *  $COPYRIGHTS$
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2011-11-11
 *
 **/
#include "common.h"
#include "drivers.h"
#if defined(HAVE_SCOTCH)
#include <scotch.h>
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Import a matrix file into a spm structure
 *
 * This function read or generate a sparse matrix from a file to store it into a
 * spm structure. The different formats accepted by this driver are described by
 * the driver field.
 *
 *******************************************************************************
 *
 * @param[in] driver
 *          This defines the driver to use to create the spm structure.
 *          = PastixDriverRSA
 *          = PastixDriverHB
 *          = PastixDriverIJV
 *          = PastixDriverMM
 *          = PastixDriverLaplacian
 *          = PastixDriverXLaplacian
 *
 * @param[in] filename
 *          The name of the file that stores the matrix (see driver)
 *
 * @param[in,out] spm
 *          On entry, an allocated sparse matrix structure.
 *          On exit, the filled sparse matrix structure with the matrix from the
 *          file.
 *
 * @param[in] pastix_comm
 *          The MPI communicator on which to distribute the sparse matrix. This
 *          is also used in case of distributed formats.
 *
 *******************************************************************************/
int
pastixReadDriver( pastix_driver_t  driver,
                  char            *filename,
                  pastix_spm_t    *spm,
                  MPI_Comm         pastix_comm )
{
    spmReadDriver( driver, filename, spm, pastix_comm );
    spmConvert( PastixCSC, spm );
    return PASTIX_SUCCESS;
}
