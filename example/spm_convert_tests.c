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
#include "expand.c"

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

int spmComp( const pastix_spm_t *spm1,
             const pastix_spm_t *spm2 )
{
    pastix_int_t *colptr1, *colptr2;
    pastix_int_t *rowptr1, *rowptr2;
    int          *valptr1, *valptr2;
    pastix_int_t  i;

    if ( spm1->fmttype != PastixCSC ) {
        fprintf(stderr, "Function made to compare only two SPM matrices\n");
        return -1;
    }

    if ((spm1->mtxtype != spm2->mtxtype) ||
        (spm1->flttype != spm2->flttype) ||
        (spm1->fmttype != spm2->fmttype) ||
        (spm1->gN      != spm2->gN     ) ||
        (spm1->n       != spm2->n      ) ||
        (spm1->gnnz    != spm2->gnnz   ) ||
        (spm1->nnz     != spm2->nnz    ) ||
        (spm1->dof     != spm2->dof    ) )
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
        pastix_int_t size = spm1->nnz * (pastix_size_of( spm1->flttype ) / sizeof(int));
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
    pastix_spm_t  spm, spm2;
    pastix_driver_t driver;
    int mtxtype, baseval;
    int ret = PASTIX_SUCCESS;
    int err = 0;

    pastix_ex_getoptions( argc, argv,
                          NULL, NULL,
                          &driver, &filename );

    spmInit(&spm);
    cscReadFromFile( driver, filename, &spm, MPI_COMM_WORLD );
    dofVar(&spm);//

    free(filename);

    printf(" -- SPM Conversion Test --\n");

    /* Allocate backup */
    memcpy( &spm2, &spm, sizeof(pastix_spm_t) );
    spm2.colptr = malloc((spm.n+1)*sizeof(pastix_int_t));
    spm2.rowptr = malloc( spm.nnz *sizeof(pastix_int_t));
    if (spm.values != NULL)
        spm2.values = malloc(spm.nnz * pastix_size_of( spm.flttype ));

    printf(" Datatype: %s\n", fltnames[spm.flttype] );
    for( baseval=0; baseval<2; baseval++ )
    {
        printf(" Baseval : %d\n", baseval );
        spmBase( &spm, baseval );

        /**
         * Backup the spm
         */
        memcpy(spm2.colptr, spm.colptr, (spm.n+1)*sizeof(pastix_int_t));
        memcpy(spm2.rowptr, spm.rowptr,  spm.nnz * sizeof(pastix_int_t));
        if (spm.values != NULL) {
            memcpy(spm2.values, spm.values, spm.nnz * pastix_size_of( spm.flttype ));
        }

        for( mtxtype=PastixGeneral; mtxtype<=PastixHermitian; mtxtype++ )
        {
            if ( (mtxtype == PastixHermitian) &&
                 ((spm.flttype != PastixComplex64) && (spm.flttype != PastixComplex32)) )
            {
                continue;
            }
            spm.mtxtype  = mtxtype;
            spm2.mtxtype = mtxtype;

            printf("   Matrix type : %s\n", mtxnames[mtxtype - PastixGeneral] );

            /**
             * Test cycle CSC -> CSR -> IJV -> CSC
             */
            printf("   -- Test Conversion CSC -> CSR: ");
            ret = spmConvert( PastixCSR, &spm );
            ret = (ret != PASTIX_SUCCESS) || (spm.fmttype != PastixCSR );
            PRINT_RES(ret);

            printf("   -- Test Conversion CSR -> IJV: ");
            ret = spmConvert( PastixIJV, &spm );
            ret = (ret != PASTIX_SUCCESS) || (spm.fmttype != PastixIJV );
            PRINT_RES(ret);

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
                ret = spmComp( &spm2, &spm );
                PRINT_RES(ret);
            }

            /**
             * Test second cycle CSC -> IJV -> CSR -> CSC
             */
            printf("   -- Test Conversion CSC -> IJV: ");
            ret = spmConvert( PastixIJV, &spm );
            ret = (ret != PASTIX_SUCCESS) || (spm.fmttype != PastixIJV );
            PRINT_RES(ret);

            printf("   -- Test Conversion IJV -> CSR: ");
            ret = spmConvert( PastixCSR, &spm );
            ret = (ret != PASTIX_SUCCESS) || (spm.fmttype != PastixCSR );
            PRINT_RES(ret);

            printf("   -- Test Conversion CSR -> CSC: ");
            ret = spmConvert( PastixCSC, &spm );
            ret = (ret != PASTIX_SUCCESS) || (spm.fmttype != PastixCSC );
            PRINT_RES(ret);

            /* Check that we came back to the initial state */
            printf("   -- Check the spm after cycle : ");
            ret = spmComp( &spm2, &spm );
            PRINT_RES(ret);
        }
        printf("\n");
    }
    spmExit( &spm  );
    spmExit( &spm2 );

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
