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

#define TEST(f) printf("-- Test %s: ", f);
#define CHECK printf("Check the CSC: "); ret = PASTIX_SUCCESS;
#define RES(a) if(a != PASTIX_SUCCESS){printf("FAILED\n"); \
               err++;} \
               else{printf("SUCCESS\n");}
#define RES_CHECK(a) if(a != PASTIX_SUCCESS){printf("FAILED\n"); \
                     err= -1; break;} \
                     else{printf("SUCCESS\n");}

int
s_compareAvals(void *avals1, void *avals2, pastix_int_t nnz)
{
    float *valptr1 = (float*)avals1;
    float *valptr2 = (float*)avals2;
    float epsilon;
    pastix_int_t i;

    /* TODO LAPACKE_dlamch_work */
    epsilon = 1E-15;

    for (i = 0; i < nnz; i++, valptr1++, valptr2++)
    {
        if ( fabs(*valptr1 - *valptr2) > epsilon )
        {
            return PASTIX_ERR_BADPARAMETER;
        }
    }

    return PASTIX_SUCCESS;
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
            return PASTIX_ERR_BADPARAMETER;
        }
    }

    return PASTIX_SUCCESS;
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
                return PASTIX_ERR_BADPARAMETER;
            }
        }
    }

    return PASTIX_SUCCESS;
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
                return PASTIX_ERR_BADPARAMETER;
            }
        }
    }

    return PASTIX_SUCCESS;
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
    int ret = PASTIX_SUCCESS;
    int err = 0;

    if( argc > 1 ) {
        filename = argv[1];
    }
    else {
        filename = "z:20:20:20";
    }

    /* Generating a 3D laplacian */
    genLaplacian( filename, &csc );

    printf("\n");
    if(csc.flttype == PastixFloat)
    {
        printf("datatype: PastixFloat\n");
    }
    else if(csc.flttype == PastixDouble)
    {
        printf("datatype: PastixDouble\n");
    }
    else if(csc.flttype == PastixComplex32)
    {
        printf("datatype: PastixComplex32\n");
    }
    else if(csc.flttype == PastixComplex64)
    {
        printf("datatype: PastixComplex64\n");
    }
    else
    {
        printf("datatype: PastixPattern\n");
    }

    printf("CSC of size %ld with %ld non-zero\n",(long)csc.gN,(long)csc.nnz);

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
        TEST("spmConvertCSC2CSR")
        ret = spmConvert( PastixCSR, &csc );
        RES(ret)

        /*-------------------------------------------------*/
        TEST("spmConvertCSR2IJV")
        ret = spmConvert( PastixIJV, &csc );
        RES(ret)

        /*-------------------------------------------------*/
        TEST("spmConvertIJV2CSC")
        ret = spmConvert( PastixCSC, &csc );
        RES(ret)

        /* intermediate check of the csc */
        /*-------------------------------------------------*/
        CHECK
        if (csc.n != n)
        {
            ret = PASTIX_ERR_BADPARAMETER;
        }
        if (csc.nnz != nnz)
        {
            ret = PASTIX_ERR_BADPARAMETER;
        }

        for (i = 0; i <= csc.n; i++)
        {
            if (csc.colptr[i] != colptr[i])
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }

        for (i = 0; i < csc.nnz; i++)
        {
            if (csc.rows[i] != rows[i])
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }

        /* check avals */
        if(csc.flttype == PastixFloat)
        {
            if (s_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixDouble)
        {
            if (d_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixComplex32)
        {
            if (c_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixComplex64)
        {
            if (z_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }
        RES_CHECK(ret)
        /* end intermediate checking of the csc */

        /*-------------------------------------------------*/
        TEST("spmConvertCSC2IJV")
        ret = spmConvert( PastixIJV, &csc );
        RES(ret)

        /*-------------------------------------------------*/
        TEST("spmConvertIJV2CSR")
        ret = spmConvert( PastixCSR, &csc );
        RES(ret)

        /*-------------------------------------------------*/
        TEST("spmConvertCSR2CSC")
        ret = spmConvert( PastixCSC, &csc );
        RES(ret)

        /* final check of the csc */
        /*-------------------------------------------------*/
        CHECK
        if (csc.n != n)
        {
            ret = PASTIX_ERR_BADPARAMETER;
        }
        if (csc.nnz != nnz)
        {
            ret = PASTIX_ERR_BADPARAMETER;
        }

        /* check colptr */
        for (i = 0; i <= csc.n; i++)
        {
            if (csc.colptr[i] != colptr[i])
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }

        /* check rows */
        for (i = 0; i < csc.nnz; i++)
        {
            if (csc.rows[i] != rows[i])
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }

        /* check avals */
        if(csc.flttype == PastixFloat)
        {
            if (s_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixDouble)
        {
            if (d_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixComplex32)
        {
            if (c_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }else if(csc.flttype == PastixComplex64)
        {
            if (z_compareAvals(csc.avals,avals,csc.nnz) > 0)
            {
                ret = PASTIX_ERR_BADPARAMETER;
            }
        }
        RES_CHECK(ret)
        /*-------------------------------------------------*/

        /* testing the rhs generator */
        /* first attempt with NULL RHS ans general csc*/
        TEST("genRHS General")
        csc.mtxtype = PastixGeneral;
        ret = genRHS(&csc,&rhs);
        RES(ret)

        /* second attempt with already filled-in RHS and hermitian csc*/
        if(csc.flttype == PastixComplex32 || csc.flttype == PastixComplex64)
        {
            csc.mtxtype=PastixHermitian;
            TEST("genRHS Hermitian")
            ret = genRHS(&csc,&rhs);
            RES(ret)
        }

        /* last try with already filled-in RHS and symetric csc*/
        csc.mtxtype=PastixSymmetric;
        TEST("genRHS Symmetric")
        ret = genRHS(&csc,&rhs);
        RES(ret)

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

    if(err==0)
    {
        printf("\n  Result: everything is OK\n\n");
        ret = PASTIX_SUCCESS;
    }
    else if(err==1)
    {
        printf("\n  Result: %d test failed\n\n",err);
        ret = PASTIX_ERR_BADPARAMETER;
    }
    else if(err==-1)
    {
        printf("\n  Abort: matrix is not correct\n\n");
        ret = PASTIX_ERR_BADPARAMETER;
    }
    else
    {
        printf("\n  Result: %d tests failed\n\n",err);
        ret = PASTIX_ERR_BADPARAMETER;
    }

    return ret;
}
