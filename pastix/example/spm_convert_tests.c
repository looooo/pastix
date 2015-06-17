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
#include <csc.h>
#include "../matrix_drivers/drivers.h"

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

int cscComp( const pastix_csc_t *csc1,
             const pastix_csc_t *csc2 )
{
    pastix_int_t *colptr1, *colptr2;
    pastix_int_t *rowptr1, *rowptr2;
    int          *valptr1, *valptr2;
    pastix_int_t  i;

    if ( csc1->fmttype != PastixCSC ) {
        fprintf(stderr, "Function made to compare only two CSC matrices\n");
        return -1;
    }

    if ((csc1->mtxtype != csc2->mtxtype) ||
        (csc1->flttype != csc2->flttype) ||
        (csc1->fmttype != csc2->fmttype) ||
        (csc1->gN      != csc2->gN     ) ||
        (csc1->n       != csc2->n      ) ||
        (csc1->gnnz    != csc2->gnnz   ) ||
        (csc1->nnz     != csc2->nnz    ) ||
        (csc1->dof     != csc2->dof    ) )
    {
        return 1;
    }

    colptr1 = csc1->colptr;
    colptr2 = csc2->colptr;
    for (i=0; i<=csc1->n; i++, colptr1++, colptr2++) {
        if (*colptr1 != *colptr2 ) {
            return 2;
        }
    }

    rowptr1 = csc1->rowptr;
    rowptr2 = csc2->rowptr;
    for (i=0; i<csc1->nnz; i++, rowptr1++, rowptr2++) {
        if (*rowptr1 != *rowptr2 ) {
            return 3;
        }
    }

    /* Check values */
    if (csc1->values != NULL) {
        pastix_int_t size = csc1->nnz * (pastix_size_of( csc1->flttype ) / sizeof(int));
        valptr1 = (int*)(csc1->values);
        valptr2 = (int*)(csc2->values);
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
    pastix_csc_t  csc, csc2;
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
    memcpy( &csc2, &csc, sizeof(pastix_csc_t) );
    csc2.colptr = malloc((csc.n+1)*sizeof(pastix_int_t));
    csc2.rowptr = malloc( csc.nnz *sizeof(pastix_int_t));
    if (csc.values != NULL)
        csc2.values = malloc(csc.nnz * pastix_size_of( csc.flttype ));

    printf(" Datatype: %s\n", fltnames[csc.flttype] );
    for( baseval=0; baseval<2; baseval++ )
    {
        printf(" Baseval : %d\n", baseval );
        spmBase( &csc, baseval );

        /**
         * Backup the csc
         */
        memcpy(csc2.colptr, csc.colptr, (csc.n+1)*sizeof(pastix_int_t));
        memcpy(csc2.rowptr, csc.rowptr,  csc.nnz * sizeof(pastix_int_t));
        if (csc.values != NULL) {
            memcpy(csc2.values, csc.values, csc.nnz * pastix_size_of( csc.flttype ));
        }

        for( mtxtype=PastixGeneral; mtxtype<=PastixHermitian; mtxtype++ )
        {
            if ( (mtxtype == PastixHermitian) &&
                 ((csc.flttype != PastixComplex64) && (csc.flttype != PastixComplex32)) )
            {
                continue;
            }
            csc.mtxtype  = mtxtype;
            csc2.mtxtype = mtxtype;

            printf("   Matrix type : %s\n", mtxnames[mtxtype - PastixGeneral] );

            /**
             * Test cycle CSC -> CSR -> IJV -> CSC
             */
            printf("   -- Test Conversion CSC -> CSR: ");
            ret = spmConvert( PastixCSR, &csc );
            ret = (ret != PASTIX_SUCCESS) || (csc.fmttype != PastixCSR );
            PRINT_RES(ret);

            printf("   -- Test Conversion CSR -> IJV: ");
            ret = spmConvert( PastixIJV, &csc );
            ret = (ret != PASTIX_SUCCESS) || (csc.fmttype != PastixIJV );
            PRINT_RES(ret);

            printf("   -- Test Conversion IJV -> CSC: ");
            ret = spmConvert( PastixCSC, &csc );
            ret = (ret != PASTIX_SUCCESS) || (csc.fmttype != PastixCSC );
            PRINT_RES(ret);

            /**
             * Check that we came back to the initial state.
             * Do not check if Symmetric or Hermitian due to transposition made
             * in the function.
             */
            if (mtxtype == PastixGeneral) {
                printf("   -- Check the csc after cycle : ");
                ret = cscComp( &csc2, &csc );
                PRINT_RES(ret);
            }

            /**
             * Test second cycle CSC -> IJV -> CSR -> CSC
             */
            printf("   -- Test Conversion CSC -> IJV: ");
            ret = spmConvert( PastixIJV, &csc );
            ret = (ret != PASTIX_SUCCESS) || (csc.fmttype != PastixIJV );
            PRINT_RES(ret);

            printf("   -- Test Conversion IJV -> CSR: ");
            ret = spmConvert( PastixCSR, &csc );
            ret = (ret != PASTIX_SUCCESS) || (csc.fmttype != PastixCSR );
            PRINT_RES(ret);

            printf("   -- Test Conversion CSR -> CSC: ");
            ret = spmConvert( PastixCSC, &csc );
            ret = (ret != PASTIX_SUCCESS) || (csc.fmttype != PastixCSC );
            PRINT_RES(ret);

            /* Check that we came back to the initial state */
            printf("   -- Check the csc after cycle : ");
            ret = cscComp( &csc2, &csc );
            PRINT_RES(ret);
        }
        printf("\n");
    }
    spmExit( &csc  );
    spmExit( &csc2 );

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
