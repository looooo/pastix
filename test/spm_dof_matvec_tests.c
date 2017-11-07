/**
 *
 * @file spm_dof_matvec_tests.c
 *
 * Tests and validate the spmMatVec routines when the spm hold constant and/or variadic dofs.
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
#include "spm.h"

int z_spm_matvec_check( int trans, const pastix_spm_t *spm );
int c_spm_matvec_check( int trans, const pastix_spm_t *spm );
int d_spm_matvec_check( int trans, const pastix_spm_t *spm );
int s_spm_matvec_check( int trans, const pastix_spm_t *spm );

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
    pastix_spm_t    original, *spm;
    pastix_driver_t driver;
    char *filename;
    int spmtype, mtxtype, fmttype, baseval;
    int ret = PASTIX_SUCCESS;
    int err = 0;
    int i, dofmax = 3;

    /**
     * Get options from command line
     */
    pastixGetOptions( argc, argv,
                      NULL, NULL,
                      NULL, &driver, &filename );

    spmReadDriver( driver, filename, &original, MPI_COMM_WORLD );
    free(filename);

    if ( original.flttype == PastixPattern ) {
        spmGenFakeValues( &original );
    }

    spmtype = original.mtxtype;
    printf(" -- SPM Matrix-Vector Test --\n");

    printf(" Datatype: %s\n", fltnames[original.flttype] );
    for( i=0; i<1; i++ )
    {
        for( mtxtype=PastixGeneral; mtxtype<=PastixHermitian; mtxtype++ )
        {
            if ( (mtxtype == PastixHermitian) &&
                 ( ((original.flttype != PastixComplex64) && (original.flttype != PastixComplex32)) ||
                   (spmtype != PastixHermitian) ) )
            {
                continue;
            }
            if ( (mtxtype != PastixGeneral) &&
                 (spmtype == PastixGeneral) )
            {
                continue;
            }
            original.mtxtype = mtxtype;

            for( baseval=0; baseval<2; baseval++ )
            {
                spmBase( &original, baseval );

                printf("   Matrix type : %s\n", mtxnames[mtxtype - PastixGeneral] );
                printf("   -- Test Matrix * Vector : ");

                /* For now only CSC is working */
                for( fmttype=0; fmttype<1; fmttype++ )
                {
                    spmConvert( fmttype, &original );
                    spm = spmDofExtend( &original, i, dofmax );

                    switch( original.flttype ){
                    case PastixComplex64:
                        ret = z_spm_matvec_check( PastixNoTrans, spm );
                        break;

                    case PastixComplex32:
                        ret = c_spm_matvec_check( PastixNoTrans, spm );
                        break;

                    case PastixFloat:
                        ret = s_spm_matvec_check( PastixNoTrans, spm );
                        break;

                    case PastixDouble:
                    default:
                        ret = d_spm_matvec_check( PastixNoTrans, spm );
                    }
                    PRINT_RES(ret);
                    spmExit( spm );
                    free(spm);
                    spm = NULL;
                }
            }
        }
    }
    spmExit( &original );

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
