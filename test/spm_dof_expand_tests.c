/**
 *
 * @file spm_dof_expand_test.c
 *
 * Tests and validate the spmNorm routines when the spm hold constant and/or variadic dofs.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 **/
#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <pastix.h>
#include "../matrix_drivers/drivers.h"
#include "spm.h"

int z_spm_norm_check( const pastix_spm_t *spm );
int c_spm_norm_check( const pastix_spm_t *spm );
int d_spm_norm_check( const pastix_spm_t *spm );
int s_spm_norm_check( const pastix_spm_t *spm );

void z_spm_print_check( char *filename, const pastix_spm_t *spm );
void c_spm_print_check( char *filename, const pastix_spm_t *spm );
void d_spm_print_check( char *filename, const pastix_spm_t *spm );
void s_spm_print_check( char *filename, const pastix_spm_t *spm );

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
    int i, rc, dofmax = 3;

    /**
     * Get options from command line
     */
    pastix_ex_getoptions( argc, argv,
                          NULL, NULL,
                          &driver, &filename );

    spmReadDriver( driver, filename, &original, MPI_COMM_WORLD );
    free(filename);

    spmtype = original.mtxtype;
    printf(" -- SPM Dof Expand Test --\n");

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
                    spm = spmDofExtend( i, dofmax, &original );

                    rc = asprintf( &filename, "%d_%s_%d_%s_%s",
                                   i, fmtnames[fmttype], baseval,
                                   mtxnames[mtxtype - PastixGeneral],
                                   fltnames[spm->flttype] );

                    printf( "-- %s --\n", filename );
                    switch( spm->flttype ){
                    case PastixComplex64:
                        z_spm_print_check( filename, spm );
                        break;

                    case PastixComplex32:
                        c_spm_print_check( filename, spm );
                        break;

                    case PastixFloat:
                        s_spm_print_check( filename, spm );
                        break;

                    case PastixDouble:
                    default:
                        d_spm_print_check( filename, spm );
                    }
                    free(filename);

                    spmExit( spm );
                    free(spm);
                    spm = NULL;
                }
            }
        }
    }
    spmExit( &original );

    (void)rc;
    return EXIT_SUCCESS;
}
