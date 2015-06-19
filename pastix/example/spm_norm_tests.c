/**
 *
 * @file spm_matvec_test.c
 *
 * Tests and validate the spm_matvec routines.
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
#include <csc.h>

int z_spm_norm_check( const pastix_csc_t *spm );
int c_spm_norm_check( const pastix_csc_t *spm );
int d_spm_norm_check( const pastix_csc_t *spm );
int s_spm_norm_check( const pastix_csc_t *spm );

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
    char *filename;
    pastix_csc_t csc;
    int mtxtype, baseval;
    int ret = PASTIX_SUCCESS;
    int err = 0;

    if( argc > 1 ) {
        filename = argv[1];
    }
    else {
        filename = "d:20:20:20";
    }

    printf(" -- SPM Matrix-Vector Test --\n");

    /* Generating a 3D laplacian */
    if ( genLaplacian( filename, &csc ) != PASTIX_SUCCESS ) {
        fprintf(stderr, "Incorrect test parameter. Accept only laplacian\n");
        return EXIT_FAILURE;
    }

    printf(" Datatype: %s\n", fltnames[csc.flttype] );
    for( baseval=0; baseval<2; baseval++ )
    {
        printf(" Baseval : %d\n", baseval );
        spmBase( &csc, baseval );

        for( mtxtype=PastixGeneral; mtxtype<=PastixHermitian; mtxtype++ )
        {
            if ( (mtxtype == PastixHermitian) &&
                 ((csc.flttype != PastixComplex64) && (csc.flttype != PastixComplex32)) )
            {
                continue;
            }

            printf("   Matrix type : %s\n", mtxnames[mtxtype - PastixGeneral] );
            printf("   -- Test Matrix * Vector : ");

            switch( csc.flttype ){
            case PastixComplex64:
                ret = z_spm_norm_check( &csc );
                break;

            case PastixComplex32:
                ret = c_spm_norm_check( &csc );
                break;

            case PastixFloat:
                ret = s_spm_norm_check( &csc );
                break;

            case PastixDouble:
            default:
                ret = d_spm_norm_check( &csc );
            }
            PRINT_RES(ret);
        }
    }
    spmExit( &csc  );

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
