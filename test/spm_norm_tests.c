/**
 *
 * @file spm_norm_tests.c
 *
 * Tests and validate the spm_norm routines.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 **/
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <pastix.h>
#include <spm.h>

int z_spm_norm_check( const pastix_spm_t *spm );
int c_spm_norm_check( const pastix_spm_t *spm );
int d_spm_norm_check( const pastix_spm_t *spm );
int s_spm_norm_check( const pastix_spm_t *spm );

#define PRINT_RES(_ret_)                        \
    if(_ret_) {                                 \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* fmtnames[] = { "CSC", "CSR", "IJV" };
char* mtxnames[] = { "General", "Symmetric", "Hermitian" };

int main (int argc, char **argv)
{
    pastix_spm_t    spm;
    pastix_driver_t driver;
    char *filename;
    int spmtype, mtxtype, fmttype, baseval;
    int ret = PASTIX_SUCCESS;
    int err = 0;

    /**
     * Get options from command line
     */
    pastixGetOptions( argc, argv,
                      NULL, NULL,
                      NULL, &driver, &filename );

    spmReadDriver( driver, filename, &spm, MPI_COMM_WORLD );
    free(filename);

    if ( spm.flttype == PastixPattern ) {
        spmGenFakeValues( &spm );
    }

    spmtype = spm.mtxtype;
    printf(" -- SPM Norms Test --\n");

    for( fmttype=0; fmttype<3; fmttype++ ) {

        spmConvert( fmttype, &spm );

        for( baseval=0; baseval<2; baseval++ )
        {
            spmBase( &spm, baseval );

            for( mtxtype=PastixGeneral; mtxtype<=PastixHermitian; mtxtype++ )
            {
                if ( (mtxtype == PastixHermitian) &&
                     ( ((spm.flttype != PastixComplex64) && (spm.flttype != PastixComplex32)) ||
                       (spmtype != PastixHermitian) ) )
                {
                    continue;
                }
                if ( (mtxtype != PastixGeneral) &&
                     (spmtype == PastixGeneral) )
                {
                    continue;
                }
                spm.mtxtype = mtxtype;

                printf(" Case: %s / %s / %d / %s\n",
                       fltnames[spm.flttype],
                       fmtnames[spm.fmttype], baseval,
                       mtxnames[mtxtype - PastixGeneral] );

                switch( spm.flttype ){
                case PastixComplex64:
                    ret = z_spm_norm_check( &spm );
                    break;

                case PastixComplex32:
                    ret = c_spm_norm_check( &spm );
                    break;

                case PastixFloat:
                    ret = s_spm_norm_check( &spm );
                    break;

                case PastixDouble:
                default:
                    ret = d_spm_norm_check( &spm );
                }
                PRINT_RES(ret);
            }
        }
    }
    spmExit( &spm  );

    if( err == 0 ) {
        printf(" -- All tests PASSED --\n");
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }
}
