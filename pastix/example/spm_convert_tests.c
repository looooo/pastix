/**
 *
 * @file z_spm_convert_test.c
 *
 * Tests and validate the spm_convert routines.
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
#include <csc.h>
#include "../matrix_drivers/drivers.h"

#define RES(_ret_)                                                  \
    if((_ret_) != PASTIX_SUCCESS) {                                 \
        printf("FAILED\n");                                         \
        err++;                                                      \
    }                                                               \
    else {                                                          \
        printf("SUCCESS\n");                                        \
    }


char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* mtxnames[] = { "General", "Symmetric", "Hermitian" };

int main (int argc, char **argv)
{
    char *filename;
    pastix_csc_t csc;
    pastix_int_t *colptr = NULL;
    pastix_int_t *rowptr = NULL;
    int          *values = NULL;
    int *valptr1, *valptr2;
    pastix_int_t i, n, nnz, size;
    int mtxtype, baseval;
    int ret = PASTIX_SUCCESS;
    int err = 0;

    if( argc > 1 ) {
        filename = argv[1];
    }
    else {
        filename = "d:20:20:20";
    }

    printf(" -- SPM Conversion Test --\n");

    /* Generating a 3D laplacian */
    if ( genLaplacian( filename, &csc ) != PASTIX_SUCCESS ) {
        fprintf(stderr, "Incorrect test parameter. Accept only laplacian\n");
        return EXIT_FAILURE;
    }

    /* Allocate backup */
    n   = csc.n;
    nnz = csc.nnz;
    colptr = malloc((n+1)*sizeof(pastix_int_t));
    rowptr = malloc(nnz * sizeof(pastix_int_t));
    if (csc.values != NULL)
        values = malloc(nnz * pastix_size_of( csc.flttype ));

    printf(" Datatype: %s\n", fltnames[csc.flttype] );

    for( baseval=0; baseval<2; baseval++ )
    {
        printf(" Baseval : %d\n", baseval );
        spmBase( &csc, baseval );

        /**
         * Backup the csc
         */
        memcpy(colptr, csc.colptr, (n+1)*sizeof(pastix_int_t));
        memcpy(rowptr, csc.rowptr, nnz * sizeof(pastix_int_t));
        if (csc.values != NULL) {
            memcpy(values, csc.values, nnz * pastix_size_of( csc.flttype ));
        }

        for( mtxtype=PastixGeneral; mtxtype<=PastixHermitian; mtxtype++ )
        {
            if ( (mtxtype == PastixHermitian) &&
                 ((csc.flttype != PastixComplex64) && (csc.flttype != PastixComplex32)) )
            {
                continue;
            }
            csc.mtxtype = mtxtype;

            printf("   Matrix type : %s\n", mtxnames[mtxtype - PastixGeneral] );

            /**
             * Test cycle CSC -> CSR -> IJV -> CSC
             */
            printf("   -- Test Conversion CSC -> CSR: ");
            ret = spmConvert( PastixCSR, &csc );
            RES(ret);

            printf("   -- Test Conversion CSR -> IJV: ");
            ret = spmConvert( PastixIJV, &csc );
            RES(ret);

            printf("   -- Test Conversion IJV -> CSC: ");
            ret = spmConvert( PastixCSC, &csc );
            RES(ret);

            /**
             * Check that we came back to the initial state.
             * Do not check if Symmetric or Hermitian due to transposition made
             * in the function.
             */
            if (mtxtype == PastixGeneral) {
                printf("   -- Check the csc after cycle : ");
                ret = PASTIX_SUCCESS;
                if (csc.n != n) {
                    ret = PASTIX_ERR_BADPARAMETER;
                }
                if (csc.nnz != nnz) {
                    ret = PASTIX_ERR_BADPARAMETER;
                }

                for (i = 0; i <= csc.n; i++) {
                    if (csc.colptr[i] != colptr[i]) {
                        ret = PASTIX_ERR_BADPARAMETER;
                        break;
                    }
                }

                for (i = 0; i < csc.nnz; i++) {
                    if (csc.rowptr[i] != rowptr[i]) {
                        ret = PASTIX_ERR_BADPARAMETER;
                        break;
                    }
                }

                /* Check values */
                if (values != NULL) {
                    size = nnz * (pastix_size_of( csc.flttype ) / sizeof(int));
                    valptr1 = values;
                    valptr2 = (int*)(csc.values);
                    for (i=0; i<size; i++, valptr1++, valptr2++) {
                        if (*valptr1 != *valptr2) {
                            ret = PASTIX_ERR_BADPARAMETER;
                            break;
                        }
                    }
                }
                RES(ret);
            }

            /**
             * Test second cycle CSC -> IJV -> CSR -> CSC
             */
            printf("   -- Test Conversion CSC -> IJV: ");
            ret = spmConvert( PastixIJV, &csc );
            RES(ret);

            printf("   -- Test Conversion IJV -> CSR: ");
            ret = spmConvert( PastixCSR, &csc );
            RES(ret);

            printf("   -- Test Conversion CSR -> CSC: ");
            ret = spmConvert( PastixCSC, &csc );
            RES(ret);

            /* Check that we came back to the initial state */
            printf("   -- Check the csc after cycle : ");
            ret = PASTIX_SUCCESS;
            if (csc.n != n) {
                ret = PASTIX_ERR_BADPARAMETER;
            }
            if (csc.nnz != nnz) {
                ret = PASTIX_ERR_BADPARAMETER;
            }

            for (i = 0; i <= csc.n; i++) {
                if (csc.colptr[i] != colptr[i]) {
                    ret = PASTIX_ERR_BADPARAMETER;
                    break;
                }
            }

            for (i = 0; i < csc.nnz; i++) {
                if (csc.rowptr[i] != rowptr[i]) {
                    ret = PASTIX_ERR_BADPARAMETER;
                    break;
                }
            }

            /* Check values */
            if (values != NULL) {
                size = nnz * (pastix_size_of( csc.flttype ) / sizeof(int));
                valptr1 = values;
                valptr2 = (int*)(csc.values);
                for (i=0; i<size; i++, valptr1++, valptr2++) {
                    if (*valptr1 != *valptr2) {
                        ret = PASTIX_ERR_BADPARAMETER;
                        break;
                    }
                }
            }
            RES(ret);
        }
        printf("\n");
    }
    spmExit( &csc );

    free( colptr ); free( rowptr );
    if ( values ) free( values );

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
