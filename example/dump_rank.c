/**
 * @file dump_rank.c
 *
 * @brief A compression example that factorizes the matrix with the Just-In-Time
 * strategy and Rank-Revealing kernels.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2023-07-21
 *
 * @ingroup pastix_examples
 * @code
 *
 */
#include "../common/common.h"
#include "../blend/solver.h"
#include "../sopalin/coeftab.h"
#include <lapacke.h>
#include <spm.h>
#include <pastix.h>

int
pastixSymbolRank( const SolverMatrix * const solvmtr,
                  FILE               * const stream )
{
    const SolverCblk *cblktnd;
    const SolverCblk *cblkptr;
    const SolverBlok *bloktnd;
    const SolverBlok *blokptr;
    int o, version = -1;

#if defined(PASTIX_SUPERNODE_STATS)
    version = 2;
#endif

    if ( stream == NULL ) {
        fprintf( stderr, "Missing output file\n" );
        return 0;
    }

    o = (fprintf (stream, "%d\n%ld\t%ld\t%ld\t%ld\n", /* Write file header */
                  version,
                  (long) solvmtr->cblknbr,
                  (long) solvmtr->bloknbr,
                  (long) solvmtr->nodenbr,
                  (long) solvmtr->baseval) == EOF);
    for (cblkptr = solvmtr->cblktab, cblktnd = cblkptr + solvmtr->cblknbr;
         (cblkptr < cblktnd) && (o == 0); cblkptr ++)
    {
#if defined(PASTIX_SUPERNODE_STATS)
        o = (fprintf (stream, "%ld\t%ld\t%ld\t%ld\n",
                      (long) cblkptr->sndeidx,
                      (long) cblkptr->fcolnum,
                      (long) cblkptr->lcolnum,
                      (long) (cblkptr->fblokptr - solvmtr->bloktab)) == EOF);
#else
        o = (fprintf (stream, "%ld\t%ld\t%ld\n",
                      (long) cblkptr->fcolnum,
                      (long) cblkptr->lcolnum,
                      (long) (cblkptr->fblokptr - solvmtr->bloktab)) == EOF);
#endif
    }
    for (blokptr = solvmtr->bloktab, bloktnd = blokptr + solvmtr->bloknbr;
         (blokptr < bloktnd) && (o == 0); blokptr ++)
    {
        pastix_int_t rkmax = 0;
        pastix_int_t rk = 0;

        rkmax = pastix_imin( cblk_colnbr( solvmtr->cblktab + blokptr->lcblknm ),
                             blok_rownbr( blokptr ) );
        rk = rkmax;

        if ( (blokptr->LRblock[0] != NULL) &&
             (blokptr->LRblock[0]->rk != -1) )
        {
            rk = blokptr->LRblock[0]->rk;
        }

        o = (fprintf (stream, "%ld\t%ld\t%ld\t%ld\t%ld\n",
                      (long) blokptr->frownum,
                      (long) blokptr->lrownum,
                      (long) blokptr->fcblknm,
                      (long) rk, (long)rkmax ) == EOF);
    }
    return (o);
}

int
pastixSymbolRank_last( const SolverMatrix * const solvmtr,
                       FILE               * const stream )
{
    const SolverBlok *bloktnd;
    const SolverBlok *blokptr;
    int o, version = -1;

    char *names[128];
    names[0] = "A11";
    names[1] = "A12";
    names[2] = "A22";

#if defined(PASTIX_SUPERNODE_STATS)
    version = 2;
#endif

    if ( stream == NULL ) {
        fprintf( stderr, "Missing output file\n" );
        return 0;
    }

    o = (fprintf (stream, "%d\n%ld\t%ld\t%ld\t%ld\n", /* Write file header */
                  version,
                  (long) solvmtr->cblknbr,
                  (long) solvmtr->bloknbr,
                  (long) solvmtr->nodenbr,
                  (long) solvmtr->baseval) == EOF);
    for (blokptr = solvmtr->bloktab, bloktnd = blokptr + solvmtr->bloknbr;
         (blokptr < bloktnd) && (o == 0); blokptr ++)
    {
        pastix_int_t    rkmax = 0;
        pastix_int_t    rk = 0;
        pastix_fixdbl_t ratio;
        rkmax = pastix_imin( cblk_colnbr( solvmtr->cblktab + blokptr->lcblknm ),
                             blok_rownbr( blokptr ) );

        rk = rkmax;

        if ( (blokptr->LRblock[0] != NULL) &&
             (blokptr->LRblock[0]->rk != -1) ){
            pastix_int_t size_fr, size_lr;
            rk = blokptr->LRblock[0]->rk;
            size_fr = cblk_colnbr( solvmtr->cblktab + blokptr->lcblknm ) * blok_rownbr( blokptr );
            size_lr = rk * ( cblk_colnbr( solvmtr->cblktab + blokptr->lcblknm ) + blok_rownbr( blokptr ) );
            ratio = (1.0*size_fr) / (1.0*size_lr);
        }
        else{
            ratio = 1.0;
        }

        SolverCblk *facing = solvmtr->cblktab + blokptr->fcblknm;
        o = (fprintf (stream, "%s %.3g\t%ld\t%ld\t%ld\t%ld\t%ld\n",
                      names[blokptr->inlast],
                      (double)ratio,
                      (long) blok_rownbr( blokptr ),
                      (long) cblk_colnbr( solvmtr->cblktab + blokptr->lcblknm ),
                      (long) rk, (long)rkmax, (long) facing->selevtx ) == EOF);
    }
    return (o);
}

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    spm_driver_t    driver;
    char           *filename = NULL;
    spmatrix_t     *spm, spm2;
    int             scatter = 0;
    int             check   = 1;
    int             rc      = 0;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );

    /**
     * Set some default low-rank parameters
     *
     * We set when to end as a trick to perform the analysis for low-rank, and
     * unset it before factorization
     */
    iparm[IPARM_COMPRESS_WHEN]       = PastixCompressWhenBegin;
    iparm[IPARM_COMPRESS_METHOD]     = PastixCompressMethodPQRCP;
    iparm[IPARM_COMPRESS_MIN_WIDTH]  = 128;
    iparm[IPARM_COMPRESS_MIN_HEIGHT] = 25;
    dparm[DPARM_COMPRESS_TOLERANCE]  = LAPACKE_dlamch_work( 'e' );

    /**
     * Get options from command line
     */
    pastixGetOptions( argc, argv,
                      iparm, dparm,
                      &check, &scatter, &driver, &filename );

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
        rc = 0;
    }

    /**
     * Generate a Fake values array if needed for the numerical part
     */
    if ( spm->flttype == SpmPattern ) {
        spmGenFakeValues( spm );
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
     * To be sure we perform the numerical factorization in Full-rank, and keep
     * the compressed flag for the compression afterward, we change the iparm
     * between the bcsc2tab and sopalin task.
     */
    //pastix_task_numfact( pastix_data, spm );
    iparm[IPARM_COMPRESS_WHEN]   = PastixCompressNever;
    iparm[IPARM_COMPRESS_METHOD] = PastixCompressMethodSVD;

    pastix_subtask_spm2bcsc( pastix_data, spm );
    pastix_subtask_bcsc2ctab( pastix_data );

    pastix_subtask_sopalin( pastix_data );

    core_get_rklimit = core_get_rklimit_test;

    {
        pastix_int_t total;
        FILE *f;

        total = coeftabCompress( pastix_data );
        fprintf(stdout, "Gain of %ld\n", (long)total );

        f = fopen("symbol_rank.sps", "w");
        pastixSymbolRank( pastix_data->solvmatr, f );
        fclose(f);

        f = fopen("symbol_rank_last.sps", "w");
        pastixSymbolRank_last( pastix_data->solvmatr, f );
        fclose(f);
    }

    spmExit( spm );
    free( spm );
    pastixFinalize( &pastix_data );

    return rc;
}

/**
 * @endcode
 */
