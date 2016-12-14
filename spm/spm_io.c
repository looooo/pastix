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
            return PASTIX_ERR_FILE;
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
            return PASTIX_ERR_FILE;
        }

        array[i  ] = (pastix_int_t)tmp1;
        array[i+1] = (pastix_int_t)tmp2;
        array[i+2] = (pastix_int_t)tmp3;
        break;
    case 2:
        if (2 != fscanf(stream, "%ld %ld", &tmp1, &tmp2)){
            errorPrint("spmLoad: Wrong input format");
            return PASTIX_ERR_FILE;
        }

        array[i  ] = (pastix_int_t)tmp1;
        array[i+1] = (pastix_int_t)tmp2;
        break;
    case 1:
        if (1 != fscanf(stream, "%ld", &tmp1)){
            errorPrint("spmLoad: Wrong input format");
            return PASTIX_ERR_FILE;
        }

        array[i  ] = (pastix_int_t)tmp1;
        break;
    }

    return PASTIX_SUCCESS;
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
            return PASTIX_ERR_FILE;
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
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (pastix_complex64_t)(tmp1 + I * tmp2);
        array[i+1] = (pastix_complex64_t)(tmp3 + I * tmp4);
        array[i+2] = (pastix_complex64_t)(tmp5 + I * tmp6);
        break;

    case 2:
        if (4 != fscanf(stream, "%lg %lg %lg %lg",
                        &tmp1, &tmp2, &tmp3, &tmp4)){
            errorPrint("spmLoad: Wrong input format");
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (pastix_complex64_t)(tmp1 + I * tmp2);
        array[i+1] = (pastix_complex64_t)(tmp3 + I * tmp4);
        break;

    case 1:
        if (2 != fscanf(stream, "%lg %lg",
                        &tmp1, &tmp2)){
            errorPrint("spmLoad: Wrong input format");
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (pastix_complex64_t)(tmp1 + I * tmp2);
        break;
    }

    return PASTIX_SUCCESS;
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
            return PASTIX_ERR_FILE;
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
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (pastix_complex32_t)(tmp1 + I * tmp2);
        array[i+1] = (pastix_complex32_t)(tmp3 + I * tmp4);
        array[i+2] = (pastix_complex32_t)(tmp5 + I * tmp6);
        break;

    case 2:
        if (4 != fscanf(stream, "%g %g %g %g",
                        &tmp1, &tmp2, &tmp3, &tmp4)){
            errorPrint("spmLoad: Wrong input format");
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (pastix_complex32_t)(tmp1 + I * tmp2);
        array[i+1] = (pastix_complex32_t)(tmp3 + I * tmp4);
        break;

    case 1:
        if (2 != fscanf(stream, "%g %g",
                        &tmp1, &tmp2)){
            errorPrint("spmLoad: Wrong input format");
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (pastix_complex32_t)(tmp1 + I * tmp2);
        break;
    }

    return PASTIX_SUCCESS;
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
            return PASTIX_ERR_FILE;
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
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (double)(tmp1);
        array[i+1] = (double)(tmp2);
        array[i+2] = (double)(tmp3);
        break;

    case 2:
        if (2 != fscanf(stream, "%lg %lg",
                        &tmp1, &tmp2)){
            errorPrint("spmLoad: Wrong input format");
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (double)(tmp1);
        array[i+1] = (double)(tmp2);
        break;

    case 1:
        if (1 != fscanf(stream, "%lg",
                        &tmp1)){
            errorPrint("spmLoad: Wrong input format");
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (double)(tmp1);
        break;
    }

    return PASTIX_SUCCESS;
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
            return PASTIX_ERR_FILE;
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
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (float)(tmp1);
        array[i+1] = (float)(tmp2);
        array[i+2] = (float)(tmp3);
        break;

    case 2:
        if (2 != fscanf(stream, "%g %g",
                        &tmp1, &tmp2)){
            errorPrint("spmLoad: Wrong input format");
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (float)(tmp1);
        array[i+1] = (float)(tmp2);
        break;

    case 1:
        if (1 != fscanf(stream, "%g", &tmp1)){
            errorPrint("spmLoad: Wrong input format");
            return PASTIX_ERR_FILE;
        }
        array[i  ] = (float)(tmp1);
        break;
    }

    return PASTIX_SUCCESS;
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
        return PASTIX_ERR_FILE;
    }

    *n   = (pastix_int_t)tmp1;
    *dof = (int)tmp2;

    /* Read the colptr array */
    *colptr = NULL;
    MALLOC_INTERN(*colptr, (*n)+1, pastix_int_t);
    assert(*colptr);

    rc = readArrayOfInteger( infile, *n+1, *colptr );
    if ( rc != PASTIX_SUCCESS )
        return rc;

    /* Read the rowptr array */
    nnz = (*colptr)[*n]-(*colptr)[0];
    *rowptr = NULL;
    MALLOC_INTERN(*rowptr, nnz, pastix_int_t);
    assert(*rowptr);

    rc = readArrayOfInteger( infile, nnz, *rowptr );
    if ( rc != PASTIX_SUCCESS )
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
            if ( rc != PASTIX_SUCCESS )
                return rc;
            break;

        case PastixComplex32:
            *values = malloc( nval * sizeof(pastix_complex32_t) );
            readArrayOfComplex32( infile, nval, *values );
            if ( rc != PASTIX_SUCCESS )
                return rc;
            break;

        case PastixDouble:
            *values = malloc( nval * sizeof(double) );
            readArrayOfDouble( infile, nval, *values );
            if ( rc != PASTIX_SUCCESS )
                return rc;
            break;

        case PastixFloat:
            *values = malloc( nval * sizeof(float) );
            readArrayOfFloat( infile, nval, *values );
            if ( rc != PASTIX_SUCCESS )
                return rc;
            break;
        }
    }

    if ( valtype != NULL ) {
        *valtype = ft;
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Load the spm structure from a file (internal format)
 *
 * Load the spm data structure from a file store in an internal format. For now
 * this function only load a simple csc matrix.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          On entry, an allocated spm data structure.
 *          On exit, the spm filled with the information read in the file
 *
 * @param[in] infile
 *          The opened file in which the spm is stored.
 *
 *******************************************************************************/
int
spmLoad( pastix_spm_t  *spm,
         FILE          *infile )
{
    pastix_int_t colsize, rowsize;
    char line[256], *test;
    int rc = PASTIX_SUCCESS;

    /**
     * Skip comments
     */
    do {
        test = fgets( line, 256, infile );
        if ( test != line ) {
            return PASTIX_ERR_FILE;
        }
    }
    while( line[0] == '#' );

    /**
     * Read header
     */
    {
        int version, mtxtype, flttype, fmttype, dof, layout;
        long gN, n, nnz, nnzexp;

        if ( 10 != sscanf( line, "%d %d %d %d %ld %ld %ld %d %ld %d\n",
                           &version, &mtxtype, &flttype, &fmttype,
                           &gN, &n, &nnz, &dof, &nnzexp, &layout ) ) {
            return PASTIX_ERR_FILE;
        }

        /* Handle only version 1 for now */
        if (version != 1) {
            return PASTIX_ERR_BADPARAMETER;
        }

        spm->mtxtype = mtxtype;
        spm->flttype = flttype;
        spm->fmttype = fmttype;
        spm->gN      = gN;
        spm->n       = n;
        spm->nnz     = nnz;
        spm->dof     = dof;
        spm->layout  = layout;

        spmUpdateFields( spm );

        assert( nnzexp == spm->nnzexp );
        assert( spm->gN == gN );
    }

    switch(spm->fmttype){
    case PastixCSC:
        colsize = spm->n + 1;
        rowsize = spm->nnz;
        break;
    case PastixCSR:
        colsize = spm->nnz;
        rowsize = spm->n + 1;
        break;
    case PastixIJV:
        colsize = spm->nnz;
        rowsize = spm->nnz;
        break;
    }

    /**
     * Read colptr
     */
    spm->colptr = malloc( colsize * sizeof(pastix_int_t) );
    rc = readArrayOfInteger( infile, colsize, spm->colptr );
    if (rc != PASTIX_SUCCESS ) {
        return rc;
    }

    /**
     * Read rowptr
     */
    spm->rowptr = malloc( rowsize * sizeof(pastix_int_t) );
    rc = readArrayOfInteger( infile, rowsize, spm->rowptr );
    if (rc != PASTIX_SUCCESS ) {
        return rc;
    }

    /**
     * Read dofs
     */
    if ( spm->dof > 0 ) {
        spm->dofs = NULL;
    }
    else {
        spm->dofs = malloc( (spm->n+1) * sizeof(pastix_int_t) );
        rc = readArrayOfInteger( infile, spm->n+1, spm->dofs );
        if (rc != PASTIX_SUCCESS ) {
            return rc;
        }
    }

    /**
     * Read dofs
     */
    if ( spm->gN == spm->n ) {
        spm->loc2glob = NULL;
    }
    else {
        spm->loc2glob = malloc( spm->n * sizeof(pastix_int_t) );
        rc = readArrayOfInteger( infile, spm->n, spm->dofs );
        if (rc != PASTIX_SUCCESS ) {
            return rc;
        }
    }

    /**
     * Read values
     */
    if (spm->flttype == PastixPattern ) {
        spm->values = NULL;
    }
    else {
        spm->values = malloc( spm->nnzexp * pastix_size_of( spm->flttype ) );
    }

    switch( spm->flttype ) {
    case PastixPattern:
        break;
    case PastixFloat:
        rc = readArrayOfFloat( infile, spm->nnzexp, spm->values );
        break;
    case PastixDouble:
        rc = readArrayOfDouble( infile, spm->nnzexp, spm->values );
        break;
    case PastixComplex32:
        rc = readArrayOfComplex32( infile, spm->nnzexp, spm->values );
        break;
    case PastixComplex64:
        rc = readArrayOfComplex64( infile, spm->nnzexp, spm->values );
        break;
    }

    return rc;
}


static inline int
writeArrayOfComplex64( FILE               *outfile,
                       pastix_int_t        n,
                       pastix_complex64_t *array )
{
    pastix_int_t i;

    /* Write 4 by 4 */
    for (i=0; i<n; i++)
    {
        fprintf(outfile, "%e %e ", creal(array[i]), cimag(array[i]));
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    return PASTIX_SUCCESS;
}

static inline int
writeArrayOfComplex32( FILE               *outfile,
                       pastix_int_t        n,
                       pastix_complex32_t *array )
{
    pastix_int_t i;

    /* Write 4 by 4 */
    for (i=0; i<n; i++)
    {
        fprintf(outfile, "%e %e ", crealf(array[i]), cimagf(array[i]));
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    return PASTIX_SUCCESS;
}

static inline int
writeArrayOfDouble( FILE         *outfile,
                    pastix_int_t  n,
                    double       *array )
{
    pastix_int_t i;

    /* Write 4 by 4 */
    for (i=0; i<n; i++)
    {
        fprintf(outfile, "%e ", array[i]);
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    return PASTIX_SUCCESS;
}

static inline int
writeArrayOfFloat( FILE         *outfile,
                   pastix_int_t  n,
                   float        *array )
{
    pastix_int_t i;

    /* Write 4 by 4 */
    for (i=0; i<n; i++)
    {
        fprintf(outfile, "%e ", array[i]);
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    return PASTIX_SUCCESS;
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

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Save the spm structure into a file (internal format)
 *
 * Save the spm data structure into a file and stored in an internal format. For
 * now this function only save a simple csc matrix.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The sparse matrix to write into the file.
 *
 * @param[in] outfile
 *          The opened file in which to store the spm.
 *
 *******************************************************************************/
int
spmSave( pastix_spm_t *spm,
         FILE         *outfile )
{
    pastix_int_t i, colsize, rowsize;

    /**
     * Write header
     */
    fprintf( outfile,
             "# version mtxtype flttype fmttype gN n nnz dof nnzexp layout\n"
             "%d %d %d %d %ld %ld %ld %d %ld %d\n",
             1, spm->mtxtype, spm->flttype, spm->fmttype,
             (long)spm->gN, (long)spm->n, (long)spm->nnz,
             (int)spm->dof, (long)spm->nnzexp, spm->layout );

    switch(spm->fmttype){
    case PastixCSC:
        colsize = spm->n + 1;
        rowsize = spm->nnz;
        break;
    case PastixCSR:
        colsize = spm->nnz;
        rowsize = spm->n + 1;
        break;
    case PastixIJV:
        colsize = spm->nnz;
        rowsize = spm->nnz;
        break;
    default:
        colsize = 0;
        rowsize = 0;
    }

    /**
     * Write colptr
     */
    for (i=0; i<colsize; i++)
    {
        fprintf(outfile, "%ld ", (long)spm->colptr[i]);
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    /**
     * Write rowptr
     */
    for (i=0; i<rowsize; i++)
    {
        fprintf(outfile, "%ld ", (long)spm->rowptr[i]);
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    /**
     * Write dofs
     */
    if ( spm->dof <= 0 ) {
        for (i=0; i<spm->n+1; i++)
        {
            fprintf(outfile, "%ld ", (long)spm->dofs[i]);
            if (i%4 == 3) fprintf(outfile, "\n");
        }
        if ((i-1)%4 !=3) fprintf(outfile, "\n");
    }

    /**
     * Write loc2glob
     */
    if ( spm->n != spm->gN ) {
        for (i=0; i<spm->n; i++)
        {
            fprintf(outfile, "%ld ", (long)spm->loc2glob[i]);
            if (i%4 == 3) fprintf(outfile, "\n");
        }
        if ((i-1)%4 !=3) fprintf(outfile, "\n");
    }

    /**
     * Write values
     */
    switch( spm->flttype ) {
    case PastixPattern:
        break;
    case PastixFloat:
        writeArrayOfFloat( outfile, spm->nnzexp, spm->values );
        break;
    case PastixDouble:
        writeArrayOfDouble( outfile, spm->nnzexp, spm->values );
        break;
    case PastixComplex32:
        writeArrayOfComplex32( outfile, spm->nnzexp, spm->values );
        break;
    case PastixComplex64:
        writeArrayOfComplex64( outfile, spm->nnzexp, spm->values );
        break;
    }

    return PASTIX_SUCCESS;
}
