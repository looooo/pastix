/**
 *
 * @file spm_test_convert.c
 *
 * Tests and validate the spm_convert routines.
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
#include <laplacian.h>

int
s_compareAvals(void *avals1, void *avals2, pastix_int_t nnz)
{
    float *valptr1 = (float*)avals1;
    float *valptr2 = (float*)avals2;
    float epsilon;
    pastix_int_t i;

    epsilon = 1E-15;

    for (i = 0; i < nnz; i++, valptr1++, valptr2++)
    {
        if ( fabs(*valptr1 - *valptr2) > epsilon )
        {
            return 1;
        }
    }

    return 0;
}

int
d_compareAvals(void *avals1, void *avals2, pastix_int_t nnz)
{
    double *valptr1 = (double*)avals1;
    double *valptr2 = (double*)avals2;
    double epsilon;
    pastix_int_t i;

    epsilon = 1E-15;

    for (i = 0; i < nnz; i++, valptr1++, valptr2++)
    {
        if ( fabs(*valptr1 - *valptr2) > epsilon )
        {
            printf("error\n");
            return 1;
        }
    }

    return 0;
}

int
c_compareAvals(void *avals1, void *avals2, pastix_int_t nnz)
{
    pastix_complex32_t *valptr1 = (pastix_complex32_t*)avals1;
    pastix_complex32_t *valptr2 = (pastix_complex32_t*)avals2;
    float epsilon;
    pastix_int_t i;

    epsilon = 1E-15;

    for (i = 0; i < nnz; i++, valptr1++, valptr2++)
    {
        if ( fabs(creal( *valptr1 ) - creal( *valptr2 ) ) > epsilon )
        {
            if ( fabs(cimag( *valptr1 ) - cimag( *valptr2 ) ) > epsilon )
            {
                printf("error\n");
                return 1;
            }
        }
    }

    return 0;
}

int
z_compareAvals(void *avals1, void *avals2, pastix_int_t nnz)
{
    pastix_complex64_t *valptr1 = (pastix_complex64_t*)avals1;
    pastix_complex64_t *valptr2 = (pastix_complex64_t*)avals2;
    double epsilon;
    pastix_int_t i;

    epsilon = 1E-15;

    for (i = 0; i < nnz; i++, valptr1++, valptr2++)
    {
        if ( fabs(creal( *valptr1 ) - creal( *valptr2 ) ) > epsilon )
        {
            if ( fabs(cimag( *valptr1 ) - cimag( *valptr2 ) ) > epsilon )
            {
                printf("error\n");
                return 1;
            }
        }
    }

    return 0;
}

int main (int argc, char **argv)
{
    char *filename;
    pastix_csc_t csc;
    pastix_int_t *colptr = NULL;
    pastix_int_t *rows   = NULL;
    void *avals          = NULL;
    void *rhs            = NULL;
    pastix_int_t i, n, nnz, baseval;
    
    if( argc > 1 ) {
        filename = argv[1];
    }
    else {
        filename = "z:20:20:20";
    }

    /* Generating a 3D laplacian */
    genLaplacian( filename, &csc );

    printf("%ld non-zero coef\n",(long)csc.nnz);

    for( baseval=1; baseval >= 0; baseval-- )
    {
        csc.mtxtype = PastixGeneral;

        printf("Basing at %d\n",(int)baseval);
        /* the laplacian generate a 1 based matrix */
        if(baseval == 0)
        {
            for(i=0;i<=csc.n;i++)
            {
                csc.colptr[i]-=1;
            }
            for(i=0;i<csc.nnz;i++)
            {
                csc.rows[i]-=1;
            }
        }

        /* Backup the csc */
        printf("Backup the csc\n");
        colptr = malloc((csc.n+1)*sizeof(pastix_int_t));
        rows   = malloc( csc.nnz *sizeof(pastix_int_t));
        memcpy(colptr,csc.colptr,(csc.n+1)*sizeof(pastix_int_t));
        memcpy(rows,csc.rows,(csc.nnz)*sizeof(pastix_int_t));
        n=csc.n;
        nnz=csc.nnz;
        if(csc.flttype == PastixFloat)
        {
            avals  = malloc(nnz       *sizeof(float));
            memcpy(avals,csc.avals,nnz*sizeof(float));
        }
        else if(csc.flttype == PastixDouble)
        {
            avals  = malloc(nnz       *sizeof(double));
            memcpy(avals,csc.avals,nnz*sizeof(double));
        }
        else if(csc.flttype == PastixComplex32)
        {
            avals  = malloc(nnz       *sizeof(pastix_complex32_t));
            memcpy(avals,csc.avals,nnz*sizeof(pastix_complex32_t));
        }
        else if(csc.flttype == PastixComplex64)
        {
            avals  = malloc(nnz       *sizeof(pastix_complex64_t));
            memcpy(avals,csc.avals,nnz*sizeof(pastix_complex64_t));
        }
        else
        {
            avals = NULL;
        }

        /* testing conversion routines */
        printf("CSC2CSR: ");
        if(spmConvert( PastixCSR, &csc ) != PASTIX_SUCCESS )
        {
            printf("error\n");
            return PASTIX_ERR_BADPARAMETER;
        }
        printf("OK\n");

        /*-------------------------------------------------*/
        printf("CSR2IJV: ");
        if(spmConvert( PastixIJV, &csc ) != PASTIX_SUCCESS )
        {
            printf("error\n");
            return PASTIX_ERR_BADPARAMETER;
        }
        printf("OK\n");

        /*-------------------------------------------------*/
        printf("IJV2CSC: ");
        if(spmConvert( PastixCSC, &csc ) != PASTIX_SUCCESS )
        {
            printf("error\n");
            return PASTIX_ERR_BADPARAMETER;
        }
        printf("OK\n");

        /* intermediate check of the csc */
        /*-------------------------------------------------*/
        printf("size: ");
        if (csc.n != n)
        {
            printf("error: n=%ld, csc.n=%ld \n",(long)n,(long)csc.n);
            return PASTIX_ERR_BADPARAMETER;
        }
        if (csc.nnz != nnz)
        {
            printf("error: nnz=%ld, csc.nnz=%ld \n",(long)nnz,(long)csc.nnz);
            return PASTIX_ERR_BADPARAMETER;
        }
        printf("   OK\n");

        printf("colptr: ");
        for (i = 0; i <= csc.n; i++)
        {
            if (csc.colptr[i] != colptr[i])
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }
        printf(" OK\n");

        printf("rows: ");
        for (i = 0; i < csc.nnz; i++)
        {
            if (csc.rows[i] != rows[i])
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }
        printf("   OK\n");

        /* check avals */
        printf("avals: ");
        if(csc.flttype == PastixFloat)
        {
            if (s_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixDouble)
        {
            if (d_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixComplex32)
        {
            if (c_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixComplex64)
        {
            if (z_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }
        printf("  OK\n");
        /* end intermediate checking of the csc */

        /*-------------------------------------------------*/
        printf("CSC2IJV: ");
        if(spmConvert( PastixIJV, &csc ) != PASTIX_SUCCESS )
        {
            printf("error\n");
            return PASTIX_ERR_BADPARAMETER;
        }
        printf("OK\n");

        /*-------------------------------------------------*/
        printf("IJV2CSR: ");
        if(spmConvert( PastixCSR, &csc ) != PASTIX_SUCCESS )
        {
            printf("error\n");
            return PASTIX_ERR_BADPARAMETER;
        }
        printf("OK\n");

        /*-------------------------------------------------*/
        printf("CSR2CSC: ");
        if(spmConvert( PastixCSC, &csc ) != PASTIX_SUCCESS )
        {
            printf("error\n");
            return PASTIX_ERR_BADPARAMETER;
        }
        printf("OK\n");

        /* final check of the csc */
        /*-------------------------------------------------*/
        printf("size: ");
        if (csc.n != n)
        {
            printf("error: n=%ld, csc.n=%ld \n",(long)n,(long)csc.n);
            return PASTIX_ERR_BADPARAMETER;
        }
        if (csc.nnz != nnz)
        {
            printf("error: nnz=%ld, csc.nnz=%ld \n",(long)nnz,(long)csc.nnz);
            return PASTIX_ERR_BADPARAMETER;
        }
        printf("   OK\n");

        /* check colptr */
        printf("colptr: ");
        for (i = 0; i <= csc.n; i++)
        {
            if (csc.colptr[i] != colptr[i])
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }
        printf(" OK\n");

        /* check rows */
        printf("rows: ");
        for (i = 0; i < csc.nnz; i++)
        {
            if (csc.rows[i] != rows[i])
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }
        printf("   OK\n");

        /* check avals */
        printf("avals: ");
        if(csc.flttype == PastixFloat)
        {
            if (s_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixDouble)
        {
            if (d_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixComplex32)
        {
            if (c_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixComplex64)
        {
            if (z_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }
        printf("  OK\n");
        /*-------------------------------------------------*/

        /* testing the rhs generator */
        printf("genRHS: ");
        /* first attempt with NULL RHS ans general csc*/
        csc.mtxtype = PastixGeneral;
        if (genRHS(&csc,&rhs) != PASTIX_SUCCESS)
        {
            printf("error\n");
            return PASTIX_ERR_BADPARAMETER;
        }

        /* second attempt with already filled-in RHS and hermitian csc*/
        if(csc.flttype == PastixComplex32 || csc.flttype == PastixComplex64)
        {
            csc.mtxtype=PastixHermitian;
            if (genRHS(&csc,&rhs) != PASTIX_SUCCESS)
            {
                printf("error\n");
                return PASTIX_ERR_BADPARAMETER;
            }
        }

        /* last try with already filled-in RHS and symetric csc*/
        csc.mtxtype=PastixSymmetric;
        if (genRHS(&csc,&rhs) != PASTIX_SUCCESS)
        {
            printf("error\n");
            return PASTIX_ERR_BADPARAMETER;
        }
        printf(" OK\n");

        /* Free memory */
        free(colptr);
        free(rows);
        free(rhs);
        rhs = NULL;
        if(avals != NULL)
            free(avals);
    /* end of the baseval loop */
    }
    
    free(csc.colptr);
    free(csc.rows);
    if(csc.avals != NULL)
        free(csc.avals);

    return 0;
}