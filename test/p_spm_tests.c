/**
 *
 * @file p_spm_tests.c
 *
 * Tests and validate the spm_convert routines.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <pastix.h>
#include <common.h>
#include <spm.h>
#include "cblas.h"
#include "lapacke.h"
#include <p_spm.h>
#include "blend/solver.h"

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
void
p_spm_print_check( char *filename, const pastix_spm_t *spm )
{
    char *file;
    FILE *f;
    int rc;

    rc = asprintf( &file, "expand_%s_sparse_cp.dat", filename );
    if ( (f = fopen( file, "w" )) == NULL ) {
        perror("p_spm_print_check:sparse_cp");
        return;
    }
    p_spmPrint( f, spm );
    fclose(f);
    free(file);

    if ( spm->dof != 1 ) {
        pastix_spm_t *espm = p_spmExpand( spm );

        rc = asprintf( &file, "expand_%s_sparse_ucp.dat", filename );
        if ( (f = fopen( file, "w" )) == NULL ) {
            perror("p_spm_print_check:sparse_ucp");
            return;
        }
        p_spmPrint( f, espm );
        fclose(f);
        free(file);

        spmExit( espm );
        free(espm);
    }

    (void)rc;
    return;
}
