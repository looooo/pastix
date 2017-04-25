/**
 *
 * @file spm_convert_test.c
 *
 * Test and validate the spmConvert routine.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s p
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
#include <spm.h>
#include "../matrix_drivers/drivers.h"

#define PRINT_RES(_ret_)                        \
    if(_ret_ == -1) {                           \
        printf("UNDEFINED\n");                  \
    }                                           \
    else if(_ret_ > 0) {                        \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* mtxnames[] = { "General", "Symmetric", "Hermitian" };

int spmComp( const pastix_spm_t *spm1,
             const pastix_spm_t *spm2 )
{
    pastix_int_t *colptr1, *colptr2;
    pastix_int_t *rowptr1, *rowptr2;
    int          *valptr1, *valptr2;
    pastix_int_t  i;

    if ( spm1->fmttype != PastixCSC ) {
        fprintf(stderr, "Function made to compare only two SPM matrices in CSC format\n");
        return -1;
    }

    if ((spm1->mtxtype != spm2->mtxtype) ||
        (spm1->flttype != spm2->flttype) ||
        (spm1->fmttype != spm2->fmttype) ||
        (spm1->gN      != spm2->gN     ) ||
        (spm1->n       != spm2->n      ) ||
        (spm1->gnnz    != spm2->gnnz   ) ||
        (spm1->nnz     != spm2->nnz    ) ||
        (spm1->dof     != spm2->dof    ) ||
        (spm1->gNexp   != spm2->gNexp  ) ||
        (spm1->nexp    != spm2->nexp   ) ||
        (spm1->gnnzexp != spm2->gnnzexp) ||
        (spm1->nnzexp  != spm2->nnzexp ) ||
        (spm1->layout  != spm2->layout ))
    {
        return 1;
    }

    colptr1 = spm1->colptr;
    colptr2 = spm2->colptr;
    for (i=0; i<=spm1->n; i++, colptr1++, colptr2++) {
        if (*colptr1 != *colptr2 ) {
            return 2;
        }
    }

    rowptr1 = spm1->rowptr;
    rowptr2 = spm2->rowptr;
    for (i=0; i<spm1->nnz; i++, rowptr1++, rowptr2++) {
        if (*rowptr1 != *rowptr2 ) {
            return 3;
        }
    }

    /* Check values */
    if (spm1->values != NULL) {
        pastix_int_t size = spm1->nnzexp * (pastix_size_of( spm1->flttype ) / sizeof(int));
        valptr1 = (int*)(spm1->values);
        valptr2 = (int*)(spm2->values);
        for (i=0; i<size; i++, valptr1++, valptr2++) {
            if (*valptr1 != *valptr2) {
                return 4;
            }
        }
    }

    return 0;
}

int main (int argc, char **argv)
{
    char *filename;
    pastix_spm_t  spm, *spm2;
    pastix_driver_t driver;
    int mtxtype, baseval;
    int ret = PASTIX_SUCCESS;
    int err = 0;
    FILE *f;
    int rc;

    pastix_ex_getoptions( argc, argv,
                          NULL, NULL,
                          &driver, &filename );

    spmReadDriver( driver, filename, &spm, MPI_COMM_WORLD );
    free(filename);

    printf(" -- SPM Conversion Test --\n");
    spmConvert(PastixCSC, &spm);

    printf(" Datatype: %s\n", fltnames[spm.flttype] );
    for( baseval=0; baseval<2; baseval++ )
    {
        printf(" Baseval : %d\n", baseval );
        spmBase( &spm, baseval );

        /**
         * Backup the spm
         */
        spm2 = spmCopy( &spm );

        for( mtxtype=PastixGeneral; mtxtype<=PastixHermitian; mtxtype++ )
        {
            if ( (mtxtype == PastixHermitian) &&
                 ((spm.flttype != PastixComplex64) && (spm.flttype != PastixComplex32)) )
            {
                continue;
            }
            spm.mtxtype  = mtxtype;
            spm2->mtxtype = mtxtype;

            printf("   Matrix type : %s\n", mtxnames[mtxtype - PastixGeneral] );

            /**
             * Test cycle CSC -> CSR -> IJV -> CSC
             */
            rc = asprintf( &filename, "convert_b%d_%s_CSC_cycle1.dat",
                           baseval, mtxnames[mtxtype - PastixGeneral] );
            f = fopen( filename, "w" );
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion CSC -> CSR: ");
            ret = spmConvert( PastixCSR, &spm );
            ret = (ret != PASTIX_SUCCESS) || (spm.fmttype != PastixCSR );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_CSR_cycle1.dat",
                           baseval, mtxnames[mtxtype - PastixGeneral] );
            f = fopen( filename, "w" );
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion CSR -> IJV: ");
            ret = spmConvert( PastixIJV, &spm );
            ret = (ret != PASTIX_SUCCESS) || (spm.fmttype != PastixIJV );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_IJV_cycle1.dat",
                           baseval, mtxnames[mtxtype - PastixGeneral] );
            f = fopen( filename, "w" );
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion IJV -> CSC: ");
            ret = spmConvert( PastixCSC, &spm );
            ret = (ret != PASTIX_SUCCESS) || (spm.fmttype != PastixCSC );
            PRINT_RES(ret);

            /**
             * Check that we came back to the initial state.
             * Do not check if Symmetric or Hermitian due to transposition made
             * in the function.
             */
            if (mtxtype == PastixGeneral) {
                printf("   -- Check the spm after cycle : ");
                ret = spmComp( spm2, &spm );
                PRINT_RES(ret);
            }

            rc = asprintf( &filename, "convert_b%d_%s_CSC_cycle2.dat",
                           baseval, mtxnames[mtxtype - PastixGeneral] );
            f = fopen( filename, "w" );
            spmPrint( &spm, f );
            fclose(f); free(filename);

            /**
             * Test second cycle CSC -> IJV -> CSR -> CSC
             */
            printf("   -- Test Conversion CSC -> IJV: ");
            ret = spmConvert( PastixIJV, &spm );
            ret = (ret != PASTIX_SUCCESS) || (spm.fmttype != PastixIJV );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_IJV_cycle2.dat",
                           baseval, mtxnames[mtxtype - PastixGeneral] );
            f = fopen( filename, "w" );
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion IJV -> CSR: ");
            ret = spmConvert( PastixCSR, &spm );
            ret = (ret != PASTIX_SUCCESS) || (spm.fmttype != PastixCSR );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_CSR_cycle2.dat",
                           baseval, mtxnames[mtxtype - PastixGeneral] );
            f = fopen( filename, "w" );
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion CSR -> CSC: ");
            ret = spmConvert( PastixCSC, &spm );
            ret = (ret != PASTIX_SUCCESS) || (spm.fmttype != PastixCSC );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_CSC_end.dat",
                           baseval, mtxnames[mtxtype - PastixGeneral] );
            f = fopen( filename, "w" );
            spmPrint( &spm, f );
            fclose(f); free(filename);

            /* Check that we came back to the initial state */
            printf("   -- Check the spm after cycle : ");
            ret = spmComp( spm2, &spm );
            PRINT_RES(ret);
        }
        printf("\n");
        spmExit( spm2 );
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

    (void)rc;
}
