/**
 * @file old/simple.c
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * This is a simple example that:
 * reads the matrix, checks if it is correct and corrects it if needed,
 * and then runs pastix in one call.
 *
 * @version 6.3.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2023-11-07
 *
 **/
#include <pastix.h>
#include <pastix/old_api.h>
#include <spm.h>

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
    pastix_float_t *b           = NULL; /* right hand side                                           */
    pastix_int_t    iparm[IPARM_SIZE]; /* integer parameters for pastix                             */
    double          dparm[DPARM_SIZE]; /* floating parameters for pastix                            */
    char           *filename;  /* Filename(s) given by user                                 */
    int             nrhs        = 1;
    spmatrix_t     *spm, spm2;
    spm_driver_t    driver;
    void           *x, *x0 = NULL;
    size_t          size;
    int             check = 1;
    int             rc    = PASTIX_SUCCESS;

#if defined(PASTIX_WITH_MPI)
    {
        int provided;
        MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided );
    }
#endif

    /*
     * Initialize parameters to default values
     */
    iparm[IPARM_MODIFY_PARAMETER] = API_NO;
    pastix( &pastix_data, MPI_COMM_WORLD,
            -1, NULL, NULL, NULL,
            NULL, NULL, NULL, 1, iparm, dparm );

    /*
     * Update options from command line, and get the matrix filename
     */
    pastixGetOptions( argc, argv,
                      iparm, dparm,
                      &check, NULL, &driver, &filename );

    /*
     * Read Matrice
     */
    spm = malloc( sizeof( spmatrix_t ) );
    rc = spmReadDriver( driver, filename, spm );
    free( filename );
    if ( rc != SPM_SUCCESS ) {
        return rc;
    }

    spmPrintInfo( spm, stdout );

    /*
     * Check Matrix format
     */
    rc = spmCheckAndCorrect( spm, &spm2 );
    if ( rc != 0 ) {
        spmExit( spm );
        *spm = spm2;
    }

    /*
     * Generate a Fake values array if needed for the numerical part
     */
    if ( spm->flttype == SpmPattern ) {
        spmGenFakeValues( spm );
    }

    iparm[IPARM_FLOAT]    = spm->flttype;
    iparm[IPARM_MTX_TYPE] = spm->mtxtype;
    iparm[IPARM_DOF_NBR]  = spm->dof;

    /**
     * Normalize A matrix (optional, but recommended for low-rank functionality)
     */
    double normA = spmNorm( SpmFrobeniusNorm, spm );
    spmScal( 1./normA, spm );

    /*
     * Generates the b and x vector such that A * x = b
     * Compute the norms of the initial vectors if checking purpose.
     */
    size = pastix_size_of( spm->flttype ) * spm->nexp;
    x = malloc( size );
    b = malloc( size );

    if ( check )
    {
        if ( check > 1 ) {
            x0 = malloc( size );
        }
        spmGenRHS( SpmRhsRndX, nrhs, spm, x0, spm->nexp, b, spm->nexp );
        memcpy( x, b, size );
    }
    else {
        spmGenRHS( SpmRhsRndB, nrhs, spm, NULL, spm->nexp, x, spm->nexp );

        /* Apply also normalization to b vector */
        spmScalMat( 1./normA, spm, nrhs, b, spm->nexp );

        /* Save b for refinement: TODO: make 2 examples w/ or w/o refinement */
        memcpy( b, x, size );
    }

    /*
     * Call pastix
     */
    iparm[IPARM_START_TASK] = API_TASK_INIT;
    iparm[IPARM_END_TASK]   = API_TASK_CLEAN;
    rc = pastix( &pastix_data, MPI_COMM_WORLD,
                 spm->nexp, spm->colptr, spm->rowptr, spm->values,
                 NULL, NULL, x, nrhs, iparm, dparm );

    if (rc != PASTIX_SUCCESS) {
        return rc;
    }

    /*
     * Check the solution
     */
    if ( check )
    {
        rc = spmCheckAxb( dparm[DPARM_EPSILON_REFINEMENT], nrhs, spm,
                          x0, spm->nexp, b, spm->nexp, x, spm->nexp );
        if ( x0 ) {
            free( x0 );
        }
    }

    spmExit( spm );
    free( spm );
    free( b );
    free( x );

#if defined(PASTIX_WITH_MPI)
    MPI_Finalize();
#endif

    return rc;
}
