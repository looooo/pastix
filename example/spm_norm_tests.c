/**
 *
 * @file spm_norm_test.c
 *
 * Tests and validate the spm_norm routines.
 *
 * @version 5.1.0
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
#include "../matrix_drivers/drivers.h"
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
char* mtxnames[] = { "General", "Symmetric", "Hermitian" };

int main (int argc, char **argv)
{
    pastix_spm_t    spm;
    pastix_driver_t driver;
    char *filename;
    int spmtype, mtxtype, baseval;
    int ret = PASTIX_SUCCESS;
    int err = 0;

    spmInit(&spm);
    /**
     * Get options from command line
     */
    pastix_ex_getoptions( argc, argv,
                          NULL, NULL,
                          &driver, &filename );

    cscReadFromFile( driver, filename, &spm, MPI_COMM_WORLD );
    dofVar(&spm);//Test dofs
    free(filename);

    spmtype = spm.mtxtype;
    printf(" -- SPM Norms Test --\n");

    printf(" Datatype: %s\n", fltnames[spm.flttype] );
    for( baseval=0; baseval<2; baseval++ )
    {
        printf(" Baseval : %d\n", baseval );
        spmBase( &spm, baseval );
        for( mtxtype=PastixGeneral; mtxtype<=PastixHermitian; mtxtype++ )
        {
            if ( (mtxtype == PastixHermitian) &&
                 ((spm.flttype != PastixComplex64) && (spm.flttype != PastixComplex32)) )
            {
                continue;
            }
            if ( (mtxtype != PastixGeneral) &&
                 (spmtype == PastixGeneral) )
            {
                continue;
            }
            spm.mtxtype = mtxtype;
            printf("   Matrix type : %s\n", mtxnames[mtxtype - PastixGeneral] );

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
