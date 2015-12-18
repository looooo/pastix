/**
 *
 * @file spm_io.c
 *
 *  PaStiX spm routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "spm.h"

static inline int
readArrayOfInteger( FILE         *stream,
                    pastix_int_t  n,
                    pastix_int_t *array )
{
    long tmp1, tmp2, tmp3, tmp4;
    pastix_int_t i;

    /* Read 4 by 4 */
    for (i=0; i<(n-3); i+=4)
    {
        if (4 != fscanf(stream, "%ld %ld %ld %ld", &tmp1, &tmp2, &tmp3, &tmp4)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }

        array[i  ] = (pastix_int_t)tmp1;
        array[i+1] = (pastix_int_t)tmp2;
        array[i+2] = (pastix_int_t)tmp3;
        array[i+3] = (pastix_int_t)tmp4;
    }

    assert( n-i < 4 );
    switch ( n - i )
    {
    case 3:
        if (3 != fscanf(stream, "%ld %ld %ld", &tmp1, &tmp2, &tmp3)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }

        array[i  ] = (pastix_int_t)tmp1;
        array[i+1] = (pastix_int_t)tmp2;
        array[i+2] = (pastix_int_t)tmp3;
        break;
    case 2:
        if (2 != fscanf(stream, "%ld %ld", &tmp1, &tmp2)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }

        array[i  ] = (pastix_int_t)tmp1;
        array[i+1] = (pastix_int_t)tmp2;
        break;
    case 1:
        if (1 != fscanf(stream, "%ld", &tmp1)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }

        array[i  ] = (pastix_int_t)tmp1;
        break;
    }

    return EXIT_SUCCESS;
}

static inline int
readArrayOfComplex64( FILE               *stream,
                      pastix_int_t        n,
                      pastix_complex64_t *array )
{
    double tmp1, tmp2, tmp3, tmp4;
    double tmp5, tmp6, tmp7, tmp8;
    pastix_int_t i;

    /* Read 4 by 4 */
    for (i=0; i<(n-3); i+=4)
    {
        if (8 != fscanf(stream, "%lg %lg %lg %lg %lg %lg %lg %lg",
                        &tmp1, &tmp2, &tmp3, &tmp4,
                        &tmp5, &tmp6, &tmp7, &tmp8)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (pastix_complex64_t)(tmp1 + I * tmp2);
        array[i+1] = (pastix_complex64_t)(tmp3 + I * tmp4);
        array[i+2] = (pastix_complex64_t)(tmp5 + I * tmp6);
        array[i+3] = (pastix_complex64_t)(tmp7 + I * tmp8);
    }

    assert( n-i < 4 );
    switch ( n - i )
    {
    case 3:
        if (6 != fscanf(stream, "%lg %lg %lg %lg %lg %lg",
                        &tmp1, &tmp2, &tmp3, &tmp4,
                        &tmp5, &tmp6)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (pastix_complex64_t)(tmp1 + I * tmp2);
        array[i+1] = (pastix_complex64_t)(tmp3 + I * tmp4);
        array[i+2] = (pastix_complex64_t)(tmp5 + I * tmp6);
        break;

    case 2:
        if (4 != fscanf(stream, "%lg %lg %lg %lg",
                        &tmp1, &tmp2, &tmp3, &tmp4)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (pastix_complex64_t)(tmp1 + I * tmp2);
        array[i+1] = (pastix_complex64_t)(tmp3 + I * tmp4);
        break;

    case 1:
        if (2 != fscanf(stream, "%lg %lg",
                        &tmp1, &tmp2)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (pastix_complex64_t)(tmp1 + I * tmp2);
        break;
    }

    return EXIT_SUCCESS;
}

static inline int
readArrayOfComplex32( FILE               *stream,
                      pastix_int_t        n,
                      pastix_complex32_t *array )
{
    float tmp1, tmp2, tmp3, tmp4;
    float tmp5, tmp6, tmp7, tmp8;
    pastix_int_t i;

    /* Read 4 by 4 */
    for (i=0; i<(n-3); i+=4)
    {
        if (8 != fscanf(stream, "%g %g %g %g %g %g %g %g",
                        &tmp1, &tmp2, &tmp3, &tmp4,
                        &tmp5, &tmp6, &tmp7, &tmp8)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (pastix_complex32_t)(tmp1 + I * tmp2);
        array[i+1] = (pastix_complex32_t)(tmp3 + I * tmp4);
        array[i+2] = (pastix_complex32_t)(tmp5 + I * tmp6);
        array[i+3] = (pastix_complex32_t)(tmp7 + I * tmp8);
    }

    assert( n-i < 4 );
    switch ( n - i )
    {
    case 3:
        if (6 != fscanf(stream, "%g %g %g %g %g %g",
                        &tmp1, &tmp2, &tmp3, &tmp4,
                        &tmp5, &tmp6)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (pastix_complex32_t)(tmp1 + I * tmp2);
        array[i+1] = (pastix_complex32_t)(tmp3 + I * tmp4);
        array[i+2] = (pastix_complex32_t)(tmp5 + I * tmp6);
        break;

    case 2:
        if (4 != fscanf(stream, "%g %g %g %g",
                        &tmp1, &tmp2, &tmp3, &tmp4)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (pastix_complex32_t)(tmp1 + I * tmp2);
        array[i+1] = (pastix_complex32_t)(tmp3 + I * tmp4);
        break;

    case 1:
        if (2 != fscanf(stream, "%g %g",
                        &tmp1, &tmp2)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (pastix_complex32_t)(tmp1 + I * tmp2);
        break;
    }

    return EXIT_SUCCESS;
}

static inline int
readArrayOfDouble( FILE         *stream,
                   pastix_int_t  n,
                   double       *array )
{
    double tmp1, tmp2, tmp3, tmp4;
    pastix_int_t i;

    /* Read 4 by 4 */
    for (i=0; i<(n-3); i+=4)
    {
        if (4 != fscanf(stream, "%lg %lg %lg %lg",
                        &tmp1, &tmp2, &tmp3, &tmp4)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (double)(tmp1);
        array[i+1] = (double)(tmp2);
        array[i+2] = (double)(tmp3);
        array[i+3] = (double)(tmp4);
    }

    assert( n-i < 4 );
    switch ( n - i )
    {
    case 3:
        if (1 != fscanf(stream, "%lg %lg %lg",
                        &tmp1, &tmp2, &tmp3)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (double)(tmp1);
        array[i+1] = (double)(tmp2);
        array[i+2] = (double)(tmp3);
        break;

    case 2:
        if (2 != fscanf(stream, "%lg %lg",
                        &tmp1, &tmp2)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (double)(tmp1);
        array[i+1] = (double)(tmp2);
        break;

    case 1:
        if (1 != fscanf(stream, "%lg",
                        &tmp1)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (double)(tmp1);
        break;
    }

    return EXIT_SUCCESS;
}

static inline int
readArrayOfFloat( FILE         *stream,
                  pastix_int_t  n,
                  float        *array )
{
    float tmp1, tmp2, tmp3, tmp4;
    pastix_int_t i;

    /* Read 4 by 4 */
    for (i=0; i<(n-3); i+=4)
    {
        if (4 != fscanf(stream, "%g %g %g %g",
                        &tmp1, &tmp2, &tmp3, &tmp4)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (float)(tmp1);
        array[i+1] = (float)(tmp2);
        array[i+2] = (float)(tmp3);
        array[i+3] = (float)(tmp4);
    }

    assert( n-i < 4 );
    switch ( n - i )
    {
    case 3:
        if (3 != fscanf(stream, "%g %g %g",
                        &tmp1, &tmp2, &tmp3)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (float)(tmp1);
        array[i+1] = (float)(tmp2);
        array[i+2] = (float)(tmp3);
        break;

    case 2:
        if (2 != fscanf(stream, "%g %g",
                        &tmp1, &tmp2)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (float)(tmp1);
        array[i+1] = (float)(tmp2);
        break;

    case 1:
        if (1 != fscanf(stream, "%g", &tmp1)){
            errorPrint("spmLoad: Wrong input format");
            return EXIT_FAILURE;
        }
        array[i  ] = (float)(tmp1);
        break;
    }

    return EXIT_SUCCESS;
}

/*
 Function: spm_load

 Load a spm from disk.

 Fill *n*, *colptr*, *rowptr*, *values* and *dof* from *infile*.

 Parameters:
 n       - number of columns
 colptr  - First cscd starting index of each column in *ja* and *a*
 rowptr    - Row of each element in first CSCD
 values  - value of each cscd in first CSCD (can be NULL)
 dof     - Number of degrees of freedom
 outfile - Output stream.

 Return:
 PASTIX_SUCCESS

 */
int
csc_load( pastix_int_t  *n,
          pastix_int_t **colptr,
          pastix_int_t **rowptr,
          int           *valtype,
          void         **values,
          int           *dof,
          FILE          *infile )
{
    int  rc, ft;
    long tmp1, tmp2;
    pastix_int_t nnz;

    /* Read the header file */
    if (3 != fscanf(infile, "%ld %ld %d\n", &tmp1, &tmp2, &ft)) {
        errorPrint("spmLoad:line 1: Wrong input");
        return EXIT_FAILURE;
    }

    *n   = (pastix_int_t)tmp1;
    *dof = (int)tmp2;

    /* Read the colptr array */
    *colptr = NULL;
    MALLOC_INTERN(*colptr, (*n)+1, pastix_int_t);
    assert(*colptr);

    rc = readArrayOfInteger( infile, *n+1, *colptr );
    if ( rc != EXIT_SUCCESS )
        return rc;

    /* Read the rowptr array */
    nnz = (*colptr)[*n]-(*colptr)[0];
    *rowptr = NULL;
    MALLOC_INTERN(*rowptr, nnz, pastix_int_t);
    assert(*rowptr);

    rc = readArrayOfInteger( infile, nnz, *rowptr );
    if ( rc != EXIT_SUCCESS )
        return rc;

    /* Read values if values is provided and if file contains */
    if (values != NULL) {
        pastix_int_t nval = nnz * (*dof) * (*dof);
        (*values) = NULL;

        switch(ft)
        {
        case PastixComplex64:
            *values = malloc( nval * sizeof(pastix_complex64_t) );
            readArrayOfComplex64( infile, nval, *values );
            if ( rc != EXIT_SUCCESS )
                return rc;
            break;

        case PastixComplex32:
            *values = malloc( nval * sizeof(pastix_complex32_t) );
            readArrayOfComplex32( infile, nval, *values );
            if ( rc != EXIT_SUCCESS )
                return rc;
            break;

        case PastixDouble:
            *values = malloc( nval * sizeof(double) );
            readArrayOfDouble( infile, nval, *values );
            if ( rc != EXIT_SUCCESS )
                return rc;
            break;

        case PastixFloat:
            *values = malloc( nval * sizeof(float) );
            readArrayOfFloat( infile, nval, *values );
            if ( rc != EXIT_SUCCESS )
                return rc;
            break;
        }
    }

    if ( valtype != NULL ) {
        *valtype = ft;
    }
    return PASTIX_SUCCESS;
}

int
spmLoad( pastix_spm_t  *spm,
         FILE          *infile )
{
    int rc, dof = 1;

    rc = csc_load( &(spm->n),
                   &(spm->colptr),
                   &(spm->rowptr),
                   &(spm->flttype),
                   &(spm->values),
                   &dof,
                   infile );

    spm->gN = spm->n;
    spm->loc2glob = NULL;

    return rc;
}


int
csc_save( pastix_int_t  n,
          pastix_int_t *colptr,
          pastix_int_t *rowptr,
          int           ft,
          void         *values,
          int           dof,
          FILE         *outfile )
{
    pastix_int_t i;

    /* Write header N Dof FloatType */
    fprintf( outfile, "%ld %ld %d\n",
             (long)n, (long)dof,
             (values == NULL) ? 0 : ft );

    /* Write colptr */
    for (i=0; i<n+1; i++)
    {
        fprintf(outfile, "%ld ", (long)colptr[i]);
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    /* Write rowptr */
    for (i=0; i<colptr[n]-1; i++)
    {
        fprintf(outfile, "%ld ", (long)rowptr[i]);
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    /* Write the values */
    if (values != NULL)
    {
        /*         for (i=0; i<(colptr[n]-1)*dof*dof; i++) */
        /*         { */
        /* #ifdef TYPE_COMPLEX */
        /*             fprintf(outfile, "%lg %lg ", (double)(creal(values[i])), (double)(cimag(values[i]))); */
        /* #else */
        /*             fprintf(outfile, "%lg ", (double)(values[i])); */
        /* #endif */
        /*             if (i%4 == 3) fprintf(outfile, "\n"); */
        /*         } */
        /*         if ((i-1)%4 !=3) fprintf(outfile, "\n"); */
    }
    return PASTIX_SUCCESS;
}

int
spmSave( pastix_spm_t *spm,
         FILE         *outfile )
{
    return csc_save( spm->n,
                     spm->colptr,
                     spm->rowptr,
                     spm->flttype,
                     spm->values,
                     1,
                     outfile );
}
