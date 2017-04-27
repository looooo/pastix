/**
 *
 * @file spm_norm_dof_test.c
 *
 * Tests and validate the spm_norm routines when the spm hold constant and/or variadic dofs.
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
    int i, dofmax = 4;

    /**
     * Get options from command line
     */
    pastix_ex_getoptions( argc, argv,
                          NULL, NULL,
                          &driver, &filename );

    spmReadDriver( driver, filename, &original, MPI_COMM_WORLD );
    free(filename);

    spmtype = original.mtxtype;
    printf(" -- SPM Norms Dof Test --\n");

    for( i=0; i<2; i++ )
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

                for( fmttype=0; fmttype<3; fmttype++ )
                {
                    spmConvert( fmttype, &original );
                    spm = spmDofExtend( &original, i, dofmax );

                    printf(" Case: %d / %s / %s / %d / %s\n",
                           i, fltnames[spm->flttype],
                           fmtnames[spm->fmttype], baseval,
                           mtxnames[mtxtype - PastixGeneral] );

                    switch( spm->flttype ){
                    case PastixComplex64:
                        ret = z_spm_norm_check( spm );
                        break;

                    case PastixComplex32:
                        ret = c_spm_norm_check( spm );
                        break;

                    case PastixFloat:
                        ret = s_spm_norm_check( spm );
                        break;

                    case PastixDouble:
                    default:
                        ret = d_spm_norm_check( spm );
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
