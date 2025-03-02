/**
 * @file example/schur.c
 *
 * @brief Schur usage example.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Matias Hastaran
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2024-07-05
 *
 * @ingroup pastix_examples
 * @code
 *
 */
#include <pastix.h>
#include <spm.h>
#include <lapacke.h>

void
schurFactorize( pastix_coeftype_t  flttype,
                pastix_factotype_t factotype,
                pastix_int_t       N,
                void              *S,
                pastix_int_t       lds,
                int              **ipiv )
{
    int info = 0;

    assert( ipiv != NULL );
    if ( factotype == PastixFactGETRF ) {
        *ipiv = malloc( N * sizeof(int) );
    }

    switch (flttype) {
    case PastixFloat:
        switch (factotype) {
        case PastixFactPOTRF:
            info = LAPACKE_spotrf_work( LAPACK_COL_MAJOR, 'L', N, S, lds );
            break;
        case PastixFactGETRF:
            info = LAPACKE_sgetrf_work( LAPACK_COL_MAJOR, N, N, S, lds, *ipiv );
            break;
        default:
            fprintf(stderr, "Factorization type not handled by Schur example\n");
        }
        break;
    case PastixComplex32:
        switch (factotype) {
        case PastixFactPOTRF:
            info = LAPACKE_cpotrf_work( LAPACK_COL_MAJOR, 'L', N, S, lds );
            break;
        case PastixFactGETRF:
            info = LAPACKE_cgetrf_work( LAPACK_COL_MAJOR, N, N, S, lds, *ipiv );
            break;
        default:
            fprintf(stderr, "Factorization type not handled by Schur example\n");
        }
        break;
    case PastixComplex64:
        switch (factotype) {
        case PastixFactPOTRF:
            info = LAPACKE_zpotrf_work( LAPACK_COL_MAJOR, 'L', N, S, lds );
            break;
        case PastixFactGETRF:
            info = LAPACKE_zgetrf_work( LAPACK_COL_MAJOR, N, N, S, lds, *ipiv );
            break;
        default:
            fprintf(stderr, "Factorization type not handled by Schur example\n");
        }
        break;
    case PastixDouble:
        switch (factotype) {
        case PastixFactPOTRF:
            info = LAPACKE_dpotrf_work( LAPACK_COL_MAJOR, 'L', N, S, lds );
            break;
        case PastixFactGETRF:
            info = LAPACKE_dgetrf_work( LAPACK_COL_MAJOR, N, N, S, lds, *ipiv );
            break;
        default:
            fprintf(stderr, "Factorization type not handled by Schur example\n");
        }
        break;
    default:
        fprintf(stderr, "Incorrect arithmetic type\n");
    }
    if (info != 0) {
        fprintf(stderr, "Error in schurFactorize with info =%d\n", info );
    }
    return;
}

void
schurSolve( pastix_coeftype_t  flttype,
            pastix_factotype_t factotype,
            pastix_int_t       Nschur,
            pastix_int_t       NRHS,
            void              *S,
            pastix_int_t       lds,
            void              *bptr,
            pastix_int_t       ldb,
            int              **ipiv )
{
    int info = 0;

    assert(ipiv != NULL);

    switch (flttype) {
    case PastixFloat:
    {
        float *b = (float *)bptr;

        switch (factotype) {
        case PastixFactPOTRF:
            info = LAPACKE_spotrs_work( LAPACK_COL_MAJOR, 'L', Nschur, NRHS, S, lds, b, ldb );
            break;
        case PastixFactGETRF:
            info = LAPACKE_sgetrs_work( LAPACK_COL_MAJOR, 'N', Nschur, NRHS, S, lds, *ipiv, b, ldb );
            break;
        default:
            fprintf(stderr, "Factorization type not handled by Schur example\n");
        }
    }
    break;
    case PastixComplex32:
    {
        pastix_complex32_t *b = (pastix_complex32_t *)bptr;

        switch (factotype) {
        case PastixFactPOTRF:
            info = LAPACKE_cpotrs_work( LAPACK_COL_MAJOR, 'L', Nschur, NRHS, S, lds, b, ldb );
            break;
        case PastixFactGETRF:
            info = LAPACKE_cgetrs_work( LAPACK_COL_MAJOR, 'N', Nschur, NRHS, S, lds, *ipiv, b, ldb );
            break;
        default:
            fprintf(stderr, "Factorization type not handled by Schur example\n");
        }
    }
    break;
    case PastixComplex64:
    {
        pastix_complex64_t *b = (pastix_complex64_t *)bptr;

        switch (factotype) {
        case PastixFactPOTRF:
            info = LAPACKE_zpotrs_work( LAPACK_COL_MAJOR, 'L', Nschur, NRHS, S, lds, b, ldb );
            break;
        case PastixFactGETRF:
            info = LAPACKE_zgetrs_work( LAPACK_COL_MAJOR, 'N', Nschur, NRHS, S, lds, *ipiv, b, ldb );
            break;
        default:
            fprintf(stderr, "Factorization type not handled by Schur example\n");
        }
    }
    break;
    case PastixDouble:
    {
        double *b = (double *)bptr;

        switch (factotype) {
        case PastixFactPOTRF:
            info = LAPACKE_dpotrs_work( LAPACK_COL_MAJOR, 'L', Nschur, NRHS, S, lds, b, ldb );
            break;
        case PastixFactGETRF:
            info = LAPACKE_dgetrs_work( LAPACK_COL_MAJOR, 'N', Nschur, NRHS, S, lds, *ipiv, b, ldb );
            break;
        default:
            fprintf(stderr, "Factorization type not handled by Schur example\n");
        }
    }
    break;
    default:
        fprintf(stderr, "Incorrect arithmetic type\n");
    }

    if (*ipiv != NULL) {
        free( *ipiv );
        *ipiv = NULL;
    }

    if (info != 0) {
        fprintf(stderr, "Error in schurSolve with info =%d\n", info );
    }

    return;
}

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    spm_driver_t    driver;
    char           *filename = NULL;
    spmatrix_t     *spm, spm2;
    void           *x, *b, *S, *x0 = NULL;
    size_t          size;
    int             scatter = 0;
    int             check   = 1;
    int             nrhs    = 1;
    int             rc      = 0;
    pastix_int_t    nschur, lds, ldb;
    int            *ipiv = NULL;
    pastix_diag_t   diag = PastixNonUnit;
    pastix_rhs_t    Xp;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );

    /**
     * Get options from command line
     */
    pastixGetOptions( argc, argv,
                      iparm, dparm,
                      &check, &scatter, &driver, &filename );


    if ( (iparm[IPARM_FACTORIZATION] == PastixFactLDLT) ||
         (iparm[IPARM_FACTORIZATION] == PastixFactLDLH) )
    {
        fprintf(stderr, "This types of factorization (LDL^t and LDL^h) are not supported by this example.\n");
        return EXIT_FAILURE;
    }

    /**
     * Startup PaStiX
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Read the sparse matrix with the driver
     */
    spm = malloc( sizeof( spmatrix_t ) );
    if ( scatter ) {
        rc = spmReadDriverDist( driver, filename, spm, MPI_COMM_WORLD );
    }
    else {
        rc = spmReadDriver( driver, filename, spm );
    }
    free( filename );
    if ( rc != SPM_SUCCESS ) {
        pastixFinalize( &pastix_data );
        return rc;
    }

    spmPrintInfo( spm, stdout );

    rc = spmCheckAndCorrect( spm, &spm2 );
    if ( rc != 0 ) {
        spmExit( spm );
        *spm = spm2;
    }

    /**
     * Generate a Fake values array if needed for the numerical part
     */
    if ( spm->flttype == SpmPattern ) {
        spmGenFakeValues( spm );
    }

    /**
     * Initialize the schur list with the first third of the unknowns
     */
    {
        nschur = spm->gN / 3;
        /* Set to a maximum to avoid memory problem with the test */
        nschur = (nschur > 5000) ? 5000 : nschur;

        if ( nschur > 0 ) {
            pastix_int_t i;
            pastix_int_t baseval = spmFindBase(spm);
            pastix_int_t *list = (pastix_int_t*)malloc(nschur * sizeof(pastix_int_t));

            for (i=0; i<nschur; i++) {
                list[i] = i+baseval;
            }
            pastixSetSchurUnknownList( pastix_data, nschur, list );
            free( list );
        }
        iparm[IPARM_SCHUR_SOLV_MODE] = PastixSolvModeInterface;
    }

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_task_analyze( pastix_data, spm );

    /**
     * Normalize A matrix (optional, but recommended for low-rank functionality)
     */
    double normA = spmNorm( SpmFrobeniusNorm, spm );
    spmScal( 1./normA, spm );

    /**
     * Perform the numerical factorization
     */
    pastix_task_numfact( pastix_data, spm );

    /**
     * Get the Schur complement back
     */
    lds = nschur;
    S = malloc( pastix_size_of( spm->flttype ) * nschur * lds );

    rc = pastixGetSchur( pastix_data, S, lds );

    if( rc == -1 ){
        spmExit( spm );
        free( spm );
        free( S );
        pastixFinalize( &pastix_data );
        exit(0);
    }

    /**
     * Factorize the Schur complement
     */
    schurFactorize( spm->flttype, iparm[IPARM_FACTORIZATION],
                    nschur, S, lds, &ipiv );

    /**
     * Generates the b and x vector such that A * x = b
     * Compute the norms of the initial vectors if checking purpose.
     */
    size = pastix_size_of( spm->flttype ) * spm->nexp * nrhs;
    x = malloc( size );
    b = malloc( size );
    ldb = spm->nexp;

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

        /* Apply also normalization to b vectors */
        spmScalMat( 1./normA, spm, nrhs, b, spm->nexp );
    }

    /**
     * Solve the linear system Ax = (P^tLUP)x = b
     */
    /* 1- Apply P to b */
    pastixRhsInit( &Xp );
    pastix_subtask_applyorder( pastix_data, PastixDirForward,
                               spm->nexp, nrhs, x, ldb, Xp );

    /* 2- Forward solve on the non Schur complement part of the system */
    if ( iparm[IPARM_FACTORIZATION] == PastixFactPOTRF ) {
        diag = PastixNonUnit;
    }
    else if( iparm[IPARM_FACTORIZATION] == PastixFactGETRF ) {
        diag = PastixUnit;
    }

    pastix_subtask_trsm( pastix_data, PastixLeft, PastixLower, PastixNoTrans, diag, Xp );

    /* 3- Solve the Schur complement part */
    {
        void *schur_x;

        size = pastix_size_of( spm->flttype ) * nschur * nrhs;
        schur_x = malloc( size );

        pastixRhsSchurGet( pastix_data, nschur, nrhs, Xp, schur_x, nschur );

        schurSolve( spm->flttype, iparm[IPARM_FACTORIZATION],
                    nschur, nrhs, S, lds, schur_x, nschur, &ipiv );

        pastixRhsSchurSet( pastix_data, nschur, nrhs, schur_x, nschur, Xp );
        free( schur_x );
    }

    /* 4- Backward solve on the non Schur complement part of the system */
    if ( iparm[IPARM_FACTORIZATION] == PastixFactPOTRF ) {
        pastix_subtask_trsm( pastix_data,
                             PastixLeft, PastixLower, PastixConjTrans, PastixNonUnit,
                             Xp );
    }
    else if( iparm[IPARM_FACTORIZATION] == PastixFactGETRF ) {
        pastix_subtask_trsm( pastix_data,
                             PastixLeft, PastixUpper, PastixNoTrans, PastixNonUnit,
                             Xp );
    }

    /* 5- Apply P^t to x */
    pastix_subtask_applyorder( pastix_data, PastixDirBackward,
                               spm->nexp, nrhs, x, ldb, Xp );
    pastixRhsFinalize( Xp );

    if ( check )
    {
        rc = spmCheckAxb( dparm[DPARM_EPSILON_REFINEMENT], nrhs, spm, x0, spm->nexp, b, spm->nexp, x, spm->nexp );

        if ( x0 ) {
            free( x0 );
        }
    }

    spmExit( spm );
    free( spm );
    free( S );
    free( x );
    free( b );
    pastixFinalize( &pastix_data );

    return rc;
}

/**
 * @endcode
 */
