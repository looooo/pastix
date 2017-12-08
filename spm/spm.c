/**
 *
 * @file spm.c
 *
 * SParse Matrix package main routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @addtogroup pastix_spm
 * @{
 **/
#include "common.h"
#include "spm.h"

#include "z_spm.h"
#include "c_spm.h"
#include "d_spm.h"
#include "s_spm.h"
#include "p_spm.h"

#include <cblas.h>

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

static int (*conversionTable[3][3][6])(pastix_spm_t*) = {
    /* From CSC */
    {{ NULL, NULL, NULL, NULL, NULL, NULL },
     { p_spmConvertCSC2CSR,
       NULL,
       s_spmConvertCSC2CSR,
       d_spmConvertCSC2CSR,
       c_spmConvertCSC2CSR,
       z_spmConvertCSC2CSR },
     { p_spmConvertCSC2IJV,
       NULL,
       s_spmConvertCSC2IJV,
       d_spmConvertCSC2IJV,
       c_spmConvertCSC2IJV,
       z_spmConvertCSC2IJV }},
    /* From CSR */
    {{ p_spmConvertCSR2CSC,
       NULL,
       s_spmConvertCSR2CSC,
       d_spmConvertCSR2CSC,
       c_spmConvertCSR2CSC,
       z_spmConvertCSR2CSC },
     { NULL, NULL, NULL, NULL, NULL, NULL },
     { p_spmConvertCSR2IJV,
       NULL,
       s_spmConvertCSR2IJV,
       d_spmConvertCSR2IJV,
       c_spmConvertCSR2IJV,
       z_spmConvertCSR2IJV }},
    /* From IJV */
    {{ p_spmConvertIJV2CSC,
       NULL,
       s_spmConvertIJV2CSC,
       d_spmConvertIJV2CSC,
       c_spmConvertIJV2CSC,
       z_spmConvertIJV2CSC },
     { p_spmConvertIJV2CSR,
       NULL,
       s_spmConvertIJV2CSR,
       d_spmConvertIJV2CSR,
       c_spmConvertIJV2CSR,
       z_spmConvertIJV2CSR },
     { NULL, NULL, NULL, NULL, NULL, NULL }}
};

#endif /* !defined(DOXYGEN_SHOULD_SKIP_THIS) */

/**
 *******************************************************************************
 *
 * @brief Init the spm structure.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to init.
 *
 *******************************************************************************/
void
spmInit( pastix_spm_t *spm )
{
    spm->mtxtype = PastixGeneral;
    spm->flttype = PastixDouble;
    spm->fmttype = PastixCSC;

    spm->gN   = 0;
    spm->n    = 0;
    spm->gnnz = 0;
    spm->nnz  = 0;

    spm->gNexp   = 0;
    spm->nexp    = 0;
    spm->gnnzexp = 0;
    spm->nnzexp  = 0;

    spm->dof      = 1;
    spm->dofs     = NULL;
    spm->layout   = PastixColMajor;

    spm->colptr   = NULL;
    spm->rowptr   = NULL;
    spm->loc2glob = NULL;
    spm->values   = NULL;
}

/**
 *******************************************************************************
 *
 * @brief Update all the computed fields based on the static values stored.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to update.
 *
 *******************************************************************************/
void
spmUpdateComputedFields( pastix_spm_t *spm )
{

    /*
     * Compute the local expended field for multi-dofs
     */
    if ( spm->dof > 0 ) {
        spm->nexp   = spm->n   * spm->dof;
        spm->nnzexp = spm->nnz * spm->dof * spm->dof;
    }
    else {
        pastix_int_t i, j, k, dofi, dofj, baseval;
        pastix_int_t *dofptr, *colptr, *rowptr;

        baseval = spmFindBase( spm );

        colptr = spm->colptr;
        rowptr = spm->rowptr;
        dofptr = spm->dofs;

        spm->nexp = dofptr[ spm->n ] - baseval;

        spm->nnzexp = 0;
        switch(spm->fmttype)
        {
        case PastixCSR:
            /* Swap pointers to call CSC */
            colptr = spm->rowptr;
            rowptr = spm->colptr;

            pastix_attr_fallthrough;

        case PastixCSC:
            for(j=0; j<spm->n; j++, colptr++) {
                dofj = dofptr[j+1] - dofptr[j];

                for(k=colptr[0]; k<colptr[1]; k++, rowptr++) {
                    i = *rowptr - baseval;
                    dofi = dofptr[i+1] - dofptr[i];

                    spm->nnzexp += dofi * dofj;
                }
            }
            break;
        case PastixIJV:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;
                dofi = dofptr[i+1] - dofptr[i];
                dofj = dofptr[j+1] - dofptr[j];

                spm->nnzexp += dofi * dofj;
            }
        }
    }

    /* TODO: add communicator */
    spm->gN      = spm->n;
    spm->gnnz    = spm->nnz;
    spm->gNexp   = spm->nexp;
    spm->gnnzexp = spm->nnzexp;
}

/**
 *******************************************************************************
 *
 * @brief Cleanup the spm structure but do not free the spm pointer.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to free.
 *
 *******************************************************************************/
void
spmExit( pastix_spm_t *spm )
{
    if(spm->colptr != NULL)
        memFree_null(spm->colptr);
    if(spm->rowptr != NULL)
        memFree_null(spm->rowptr);
    if(spm->loc2glob != NULL)
        memFree_null(spm->loc2glob);
    if(spm->values != NULL)
        memFree_null(spm->values);
    if(spm->dofs != NULL)
        memFree_null(spm->dofs);
}

/**
 *******************************************************************************
 *
 * @brief Rebase the arrays of the spm to the given value.
 *
 * If the value is equal to the original base, then nothing is performed.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to rebase.
 *
 * @param[in] baseval
 *          The base value to use in the graph (0 or 1).
 *
 *******************************************************************************/
void
spmBase( pastix_spm_t *spm,
         int           baseval )
{
    pastix_int_t baseadj;
    pastix_int_t i, n, nnz;

    /* Parameter checks */
    if ( spm == NULL ) {
        pastix_print_error("spmBase: spm pointer is NULL");
        return;
    }
    if ( (spm->colptr == NULL) ||
         (spm->rowptr == NULL) )
    {
        pastix_print_error("spmBase: spm pointer is not correctly initialized");
        return;
    }
    if ( (baseval != 0) &&
         (baseval != 1) )
    {
        pastix_print_error("spmBase: baseval is incorrect, must be 0 or 1");
        return;
    }

    baseadj = baseval - spmFindBase( spm );
    if (baseadj == 0)
	return;

    n   = spm->n;
    nnz = spm->nnz;

    switch(spm->fmttype)
    {
    case PastixCSC:
        assert( nnz == (spm->colptr[n] - spm->colptr[0]) );

        for (i = 0; i <= n; i++) {
            spm->colptr[i] += baseadj;
        }
        for (i = 0; i < nnz; i++) {
            spm->rowptr[i] += baseadj;
        }
        break;

    case PastixCSR:
        assert( nnz == (spm->rowptr[n] - spm->rowptr[0]) );
        for (i = 0; i <= n; i++) {
            spm->rowptr[i] += baseadj;
        }
        for (i = 0; i < nnz; i++) {
            spm->colptr[i] += baseadj;
        }
        break;
    case PastixIJV:
        for (i = 0; i < nnz; i++) {
            spm->rowptr[i] += baseadj;
            spm->colptr[i] += baseadj;
        }
    }

    if (spm->loc2glob != NULL) {
        for (i = 0; i < n; i++) {
            spm->loc2glob[i] += baseadj;
        }
    }
    if (spm->dofs != NULL) {
        for (i = 0; i <= n; i++) {
            spm->dofs[i] += baseadj;
        }
    }
    return;
}

/**
 *******************************************************************************
 *
 * @brief Search the base used in the spm structure.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 ********************************************************************************
 *
 * @return  The baseval used in the given sparse matrix structure.
 *
 *******************************************************************************/
pastix_int_t
spmFindBase( const pastix_spm_t *spm )
{

    pastix_int_t i, *tmp, baseval;

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    baseval = pastix_imin( *(spm->colptr), *(spm->rowptr) );
    /*
     * if not:
     */
    if ( ( baseval != 0 ) &&
         ( baseval != 1 ) )
    {
        baseval = spm->n;
        tmp = spm->colptr;
        for(i=0; i<spm->nnz; i++, tmp++){
            baseval = pastix_imin( *tmp, baseval );
        }
    }

    return baseval;
}

/**
 *******************************************************************************
 *
 * @brief  Convert the storage format of the spm.
 *
  *******************************************************************************
 *
 * @param[in] ofmttype
 *          The output format of the sparse matrix. It must be:
 *          - PastixCSC
 *          - PastixCSR
 *          - PastixIJV
 *
 * @param[inout] spm
 *          The sparse matrix structure to convert.
 *
 ********************************************************************************
 *
 * @retval PASTIX_SUCCESS if the conversion happened successfully.
 * @retval PASTIX_ERR_BADPARAMETER if one the parameter is incorrect.
 * @retval PASTIX_ERR_NOTIMPLEMENTED if the case is not yet implemented.
 *
 *******************************************************************************/
int
spmConvert( int ofmttype, pastix_spm_t *spm )
{
    if ( conversionTable[spm->fmttype][ofmttype][spm->flttype] ) {
        if ( spm->dof != 1 ) {
            //pastix_print_error( "spmConvert: Conversion of non unique dof not yet implemented\n");
            return PASTIX_ERR_NOTIMPLEMENTED;
        }
        return conversionTable[spm->fmttype][ofmttype][spm->flttype]( spm );
    }
    else {
        return PASTIX_SUCCESS;
    }
}

/**
 *******************************************************************************
 *
 * @brief Convert the spm matrix into a dense matrix for test purpose.
 *
 * @remark DO NOT USE with large matrices.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix structure to convert.
 *
 ********************************************************************************
 *
 * @return
 *        The pointer to the allocated array storing the dense version of the
 *        matrix.
 *
 *******************************************************************************/
void *
spm2Dense( const pastix_spm_t *spm )
{
    switch (spm->flttype) {
    case PastixFloat:
        return s_spm2dense( spm );
    case PastixComplex32:
        return c_spm2dense( spm );
    case PastixComplex64:
        return z_spm2dense( spm );
    case PastixDouble:
        return d_spm2dense( spm );
    default:
        return NULL;
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the norm of the spm.
 *
 * Return the ntype norm of the sparse matrix spm.
 *
 *     spmNorm = ( max(abs(spm(i,j))), NORM = PastixMaxNorm
 *               (
 *               ( norm1(spm),         NORM = PastixOneNorm
 *               (
 *               ( normI(spm),         NORM = PastixInfNorm
 *               (
 *               ( normF(spm),         NORM = PastixFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(spm(i,j))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          - PastixMaxNorm
 *          - PastixOneNorm
 *          - PastixInfNorm
 *          - PastixFrobeniusNorm
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 ********************************************************************************
 *
 * @retval norm The norm described above. Note that for simplicity, even if the
 *              norm of single real or single complex matrix is computed with
 *              single precision, the returned norm is stored in double
 *              precision.
 * @retval -1   If the floating point of the sparse matrix is undefined.
 *
 *******************************************************************************/
double
spmNorm( pastix_normtype_t   ntype,
         const pastix_spm_t *spm )
{
    pastix_spm_t *spmtmp = (pastix_spm_t*)spm;
    double norm = -1.;

    if ( spm->dof != 1 ) {
        fprintf(stderr, "WARNING: spm expanded due to non implemented norm for non-expanded spm\n");
        spmtmp = spmExpand( spm );
    }
    switch (spm->flttype) {
    case PastixFloat:
        norm = (double)s_spmNorm( ntype, spmtmp );
        break;

    case PastixDouble:
        norm = d_spmNorm( ntype, spmtmp );
        break;

    case PastixComplex32:
        norm = (double)c_spmNorm( ntype, spmtmp );
        break;

    case PastixComplex64:
        norm = z_spmNorm( ntype, spmtmp );
        break;

    case PastixPattern:
    default:
        ;
    }

    if ( spmtmp != spm ) {
        spmExit( spmtmp );
        free(spmtmp);
    }
    return norm;
}

/**
 *******************************************************************************
 *
 * @brief Sort the subarray of edges of each vertex in a CSC or CSR format.
 *
 * Nothing is performed if IJV format is used.
 *
 * @warning This function should NOT be called if dof is greater than 1.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with subarrays of edges sorted by
 *          ascending order.
 *
 ********************************************************************************
 *
 * @retval PASTIX_SUCCESS if the sort was called
 * @retval PASTIX_ERR_BADPARAMETER if the given spm was incorrect.
 *
 *******************************************************************************/
int
spmSort( pastix_spm_t *spm )
{
    if ( spm->dof != 1 ) {
        fprintf(stderr, "WARNING: spm expanded due to non implemented sort for non-expanded spm\n");
        spm = spmExpand( spm );
    }
    switch (spm->flttype) {
    case PastixPattern:
        p_spmSort( spm );
        break;
    case PastixFloat:
        s_spmSort( spm );
        break;
    case PastixDouble:
        d_spmSort( spm );
        break;
    case PastixComplex32:
        c_spmSort( spm );
        break;
    case PastixComplex64:
        z_spmSort( spm );
        break;
    default:
        return PASTIX_ERR_BADPARAMETER;
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Merge multiple entries in a spm by summing their values together.
 *
 * The sparse matrix needs to be sorted first (see spmSort()).
 *
 * @warning Not implemented for CSR and IJV format.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the reduced sparse matrix of multiple entries were present
 *          in it. The multiple values for a same vertex are sum up together.
 *
 ********************************************************************************
 *
 * @retval >=0 the number of vertices that were merged,
 * @retval PASTIX_ERR_BADPARAMETER if the given spm was incorrect.
 *
 *******************************************************************************/
pastix_int_t
spmMergeDuplicate( pastix_spm_t *spm )
{
    if ( spm->dof != 1 ) {
        fprintf(stderr, "WARNING: spm expanded due to non implemented merge for non-expanded spm\n");
        spm = spmExpand( spm );
    }
    switch (spm->flttype) {
    case PastixPattern:
        return p_spmMergeDuplicate( spm );

    case PastixFloat:
        return s_spmMergeDuplicate( spm );

    case PastixDouble:
        return d_spmMergeDuplicate( spm );

    case PastixComplex32:
        return c_spmMergeDuplicate( spm );

    case PastixComplex64:
        return z_spmMergeDuplicate( spm );

    default:
        return PASTIX_ERR_BADPARAMETER;
    }
}

/**
 *******************************************************************************
 *
 * @brief Symmetrize the pattern of the spm.
 *
 * This routine corrects the sparse matrix structure if it's pattern is not
 * symmetric. It returns the new symmetric pattern with zeroes on the new
 * entries.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with extra entries that makes it
 *          pattern symmetric.
 *
 ********************************************************************************
 *
 * @retval >=0 the number of entries added to the matrix,
 * @retval PASTIX_ERR_BADPARAMETER if the given spm was incorrect.
 *
 *******************************************************************************/
pastix_int_t
spmSymmetrize( pastix_spm_t *spm )
{
    if ( spm->dof != 1 ) {
        fprintf(stderr, "WARNING: spm expanded due to non implemented symmetrize for non-expanded spm\n");
        spm = spmExpand( spm );
    }
    switch (spm->flttype) {
    case PastixPattern:
        return p_spmSymmetrize( spm );

    case PastixFloat:
        return s_spmSymmetrize( spm );

    case PastixDouble:
        return d_spmSymmetrize( spm );

    case PastixComplex32:
        return c_spmSymmetrize( spm );

    case PastixComplex64:
        return z_spmSymmetrize( spm );

    default:
        return PASTIX_ERR_BADPARAMETER;
    }
}

/**
 *******************************************************************************
 *
 * @brief Check the correctness of a spm.
 *
 * This routine initializes the sparse matrix to fit the PaStiX requirements. If
 * needed, the format is changed to CSC, the duplicated vertices are merged
 * together by summing their values; the graph is made symmetric for matrices
 * with unsymmetric pattern, new values are set to 0. Only the lower part is
 * kept for symmetric matrices.
 *
  *******************************************************************************
 *
 * @param[inout] spm
 *          The pointer to the sparse matrix structure to check, and correct.
 *
 *******************************************************************************
 *
 * @return if no changes have been made, the initial sparse matrix is returned,
 * otherwise a copy of the sparse matrix, fixed to meet the PaStiX requirements,
 * is returned.
 *
 *******************************************************************************/
pastix_spm_t *
spmCheckAndCorrect( pastix_spm_t *spm )
{
    pastix_spm_t *newspm = NULL;
    pastix_int_t count;

    /* Let's work on a copy */
    newspm = spmCopy( spm );

    /* PaStiX works on CSC matrices */
    spmConvert( PastixCSC, newspm );

    if ( newspm->dof != 1 ) {
        fprintf(stderr, "WARNING: newspm expanded due to missing check functions implementations\n");
        newspm = spmExpand( newspm );
    }

    /* Sort the rowptr for each column */
    spmSort( newspm );

    /* Merge the duplicated entries by summing the values */
    count = spmMergeDuplicate( newspm );
    if ( count > 0 ) {
        fprintf(stderr, "spmCheckAndCorrect: %ld entries have been merged\n", (long)count );
    }

    /*
     * If the matrix is symmetric or hermitian, we keep only the upper or lower
     * part, otherwise, we symmetrize the graph to get A+A^t, new values are set
     * to 0.
     */
    if ( newspm->mtxtype == PastixGeneral ) {
        count = spmSymmetrize( newspm );
        if ( count > 0 ) {
            fprintf(stderr, "spmCheckAndCorrect: %ld entries have been added for symmetry\n", (long)count );
        }
    }
    else {
        //spmToLower( newspm );
    }

    spmUpdateComputedFields( newspm );

    /*
     * Check if we return the new one, or the original one because no changes
     * have been made
     */
    if (( spm->fmttype != newspm->fmttype ) ||
        ( spm->nnzexp  != newspm->nnzexp  ) )
    {
        return newspm;
    }
    else {
        spmExit( newspm );
        free(newspm);
        return (pastix_spm_t*)spm;
    }
}

/**
 *******************************************************************************
 *
 * @brief Create a copy of the spm.
 *
 * Duplicate the spm data structure given as parameter. All new arrays are
 * allocated and copied from the original matrix. Both matrices need to be
 * freed.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to copy.
 *
 *******************************************************************************
 *
 * @return
 *          The copy of the sparse matrix.
 *
 *******************************************************************************/
pastix_spm_t *
spmCopy( const pastix_spm_t *spm )
{
    pastix_spm_t *newspm = (pastix_spm_t*)malloc(sizeof(pastix_spm_t));
    pastix_int_t colsize, rowsize, valsize, dofsize;

    memcpy( newspm, spm, sizeof(pastix_spm_t));

    switch(spm->fmttype){
    case PastixCSC:
        colsize = spm->n + 1;
        rowsize = spm->nnz;
        valsize = spm->nnzexp;
        dofsize = spm->n + 1;
        break;
    case PastixCSR:
        colsize = spm->nnz;
        rowsize = spm->n + 1;
        valsize = spm->nnzexp;
        dofsize = spm->n + 1;
        break;
    case PastixIJV:
    default:
        colsize = spm->nnz;
        rowsize = spm->nnz;
        valsize = spm->nnzexp;
        dofsize = spm->n + 1;
    }

    if(spm->colptr != NULL) {
        newspm->colptr = (pastix_int_t*)malloc( colsize * sizeof(pastix_int_t));
        memcpy( newspm->colptr, spm->colptr, colsize * sizeof(pastix_int_t));
    }
    if(spm->rowptr != NULL) {
        newspm->rowptr = (pastix_int_t*)malloc(rowsize * sizeof(pastix_int_t));
        memcpy( newspm->rowptr, spm->rowptr, rowsize * sizeof(pastix_int_t));
    }
    if(spm->loc2glob != NULL) {
        newspm->loc2glob = (pastix_int_t*)malloc(dofsize * sizeof(pastix_int_t));
        memcpy( newspm->loc2glob, spm->loc2glob, dofsize * sizeof(pastix_int_t));
    }
    if(spm->dofs != NULL) {
        newspm->dofs = (pastix_int_t*)malloc(dofsize * sizeof(pastix_int_t));
        memcpy( newspm->dofs, spm->dofs, dofsize * sizeof(pastix_int_t) );
    }
    if(spm->values != NULL) {
        valsize = valsize * pastix_size_of( spm->flttype );
        newspm->values = malloc(valsize);
        memcpy( newspm->values, spm->values, valsize );
    }
    return newspm;
}

/**
 *******************************************************************************
 *
 * @brief Print basic informations about the spm matrix into a given stream.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to print.
 *
 * @param[inout] stream
 *          Stream to print the spm matrix. stdout is used if stream == NULL.
 *
 *******************************************************************************/
void
spmPrintInfo( const pastix_spm_t* spm, FILE *stream )
{
    char *mtxtypestr[4] = { "General", "Symmetric", "Hermitian", "Incorrect" };
    char *flttypestr[7] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64", "Incorrect" };
    char *fmttypestr[4] = { "CSC", "CSR", "IJV", "Incorrect" };
    int mtxtype = spm->mtxtype - PastixGeneral;
    int flttype = spm->flttype - PastixPattern;
    int fmttype = spm->fmttype - PastixCSC;

    if (stream == NULL) {
        stream = stdout;
    }

    mtxtype = (mtxtype > 2 || mtxtype < 0) ? 3 : mtxtype;
    flttype = (flttype > 5 || flttype < 0) ? 6 : flttype;
    fmttype = (fmttype > 2 || fmttype < 0) ? 3 : fmttype;

    fprintf(stream,
            "  Matrix type:  %s\n"
            "  Arithmetic:   %s\n"
            "  Format:       %s\n"
            "  N:            %ld\n"
            "  nnz:          %ld\n",
            mtxtypestr[mtxtype],
            flttypestr[flttype],
            fmttypestr[fmttype],
            (long)spm->gN, (long)spm->gnnz );

    if ( spm->dof != 1 ) {
        if ( spm->dof > 1 ) {
            fprintf(stream,
                    "  Dof:          %ld\n",
                    (long)spm->dof );
        }
        else {
            fprintf(stream,
                    "  Dof:          Variadic\n" );
        }

        fprintf(stream,
                "  N expanded:   %ld\n"
                "  NNZ expanded: %ld\n",
                (long)spm->gNexp, (long)spm->gnnzexp );
    }
}

/**
 *******************************************************************************
 *
 * @brief Print an spm matrix into into a given file.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to print.
 *
 * @param[in] stream
 *          File to print the spm matrix. stdout, if stream == NULL.
 *
 *******************************************************************************/
void
spmPrint( const pastix_spm_t *spm,
          FILE               *stream )
{
    if (stream == NULL) {
        stream = stdout;
    }

    switch(spm->flttype)
    {
    case PastixPattern:
        //return p_f, spmPrint(f, spm);
        break;
    case PastixFloat:
        s_spmPrint(stream, spm);
        break;
    case PastixComplex32:
        c_spmPrint(stream, spm);
        break;
    case PastixComplex64:
        z_spmPrint(stream, spm);
        break;
    case PastixDouble:
    default:
        d_spmPrint(stream, spm);
    }
}

/**
 *******************************************************************************
 *
 * @brief Expand a multi-dof spm matrix into an spm with constant dof set to 1.
 *
 * Duplicate the spm data structure given as parameter. All new arrays are
 * allocated and copied from the original matrix. Both matrices need to be
 * freed.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to copy.
 *
 *******************************************************************************
 *
 * @return
 *          The copy of the sparse matrix.
 *
 *******************************************************************************/
pastix_spm_t *
spmExpand( const pastix_spm_t* spm )
{
    switch(spm->flttype)
    {
    case PastixPattern:
        return p_spmExpand(spm);
        break;
    case PastixFloat:
        return s_spmExpand(spm);
        break;
    case PastixComplex32:
        return c_spmExpand(spm);
        break;
    case PastixComplex64:
        return z_spmExpand(spm);
        break;
    case PastixDouble:
    default:
        return d_spmExpand(spm);
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute a matrix-vector product.
 *
 *    y = alpha * op(A) * x + beta * y, where op(A) is one of
 *
 *    op( A ) = A  or op( A ) = A' or op( A ) = conjg( A' )
 *
 *  alpha and beta are scalars, and x and y are vectors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix spm is transposed, not transposed or conjugate transposed:
 *          - PastixTrans
 *          - PastixNoTrans
 *          - PastixConjTrans
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] spm
 *          The PastixGeneral spm.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta.
 *
 * @param[inout] y
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the y vector has been computed successfully,
 * @retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmMatVec(       pastix_trans_t trans,
           const void          *alpha,
           const pastix_spm_t  *spm,
           const void          *x,
           const void          *beta,
                 void          *y )
{
    pastix_spm_t *espm = (pastix_spm_t*)spm;
    int rc = PASTIX_SUCCESS;

    if ( spm->fmttype != PastixCSC ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( spm->dof != 1 ) {
        espm = spmExpand( spm );
    }
    switch (spm->flttype) {
    case PastixFloat:
        rc = s_spmCSCMatVec( trans, alpha, espm, x, beta, y );
        break;
    case PastixComplex32:
        rc = c_spmCSCMatVec( trans, alpha, espm, x, beta, y );
        break;
    case PastixComplex64:
        rc = z_spmCSCMatVec( trans, alpha, espm, x, beta, y );
        break;
    case PastixDouble:
    default:
        rc = d_spmCSCMatVec( trans, alpha, espm, x, beta, y );
    }

    if ( spm != espm ) {
        spmExit( espm );
        free(espm);
    }
    return rc;
}

/**
 *******************************************************************************
 *
 * @brief Generate right hand side vectors associated to a given matrix.
 *
 * The vectors can be initialized randomly or to get a specific solution.
 *
 *******************************************************************************
 *
 * @param[in] type
 *          Defines how to compute the vector b.
 *          - PastixRhsOne:  b is computed such that x = 1 [ + I ]
 *          - PastixRhsI:    b is computed such that x = i [ + i * I ]
 *          - PastixRhsRndX: b is computed by matrix-vector product, such that
 *            is a random vector in the range [-0.5, 0.5]
 *          - PastixRhsRndB: b is computed randomly and x is not computed.
 *
 * @param[in] nrhs
 *          Defines the number of right hand side that must be generated.
 *
 * @param[in] spm
 *          The sparse matrix uses to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[out] x
 *          On exit, if x != NULL, then the x vector(s) generated to compute b
 *          is returned. Must be of size at least ldx * spm->n.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= spm->n.
 *
 * @param[inout] b
 *          b must be an allocated matrix of size at least ldb * nrhs.
 *          On exit, b is initialized as defined by the type parameter.
 *
 * @param[in] ldb
 *          Defines the leading dimension of b when multiple right hand sides
 *          are available. ldb >= spm->n.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the b vector has been computed successfully,
 * @retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmGenRHS( pastix_rhstype_t type, pastix_int_t nrhs,
           const pastix_spm_t  *spm,
           void                *x, pastix_int_t ldx,
           void                *b, pastix_int_t ldb )
{
    static int (*ptrfunc[4])(pastix_rhstype_t, int,
                             const pastix_spm_t *,
                             void *, int, void *, int) =
        {
            s_spmGenRHS, d_spmGenRHS, c_spmGenRHS, z_spmGenRHS
        };

    int id = spm->flttype - PastixFloat;
    if ( (id < 0) || (id > 3) ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    else {
        return ptrfunc[id](type, nrhs, spm, x, ldx, b, ldb );
    }
}

/**
 *******************************************************************************
 *
 * @brief Check the backward error, and the forward error if x0 is provided.
 *
 *******************************************************************************
 *
 * @param[in] nrhs
 *          Defines the number of right hand side that must be generated.
 *
 * @param[in] spm
 *          The sparse matrix used to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[inout] x0
 *          If x0 != NULL, the forward error is computed.
 *          On exit, x0 stores x0-x
 *
 * @param[in] ldx0
 *          Defines the leading dimension of x0 when multiple right hand sides
 *          are available. ldx0 >= spm->n.
 *
 * @param[inout] b
 *          b is a matrix of size at least ldb * nrhs.
 *          On exit, b stores Ax-b.
 *
 * @param[in] ldb
 *          Defines the leading dimension of b when multiple right hand sides
 *          are available. ldb >= spm->n.
 *
 * @param[in] x
 *          Contains the solution computed by the solver.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= spm->n.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the b vector has been computed successfully,
 * @retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmCheckAxb( pastix_int_t nrhs,
             const pastix_spm_t  *spm,
                   void *x0, pastix_int_t ldx0,
                   void *b,  pastix_int_t ldb,
             const void *x,  pastix_int_t ldx )
{
    static int (*ptrfunc[4])(int, const pastix_spm_t *,
                             void *, int, void *, int, const void *, int) =
        {
            s_spmCheckAxb, d_spmCheckAxb, c_spmCheckAxb, z_spmCheckAxb
        };

    int id = spm->flttype - PastixFloat;
    if ( (id < 0) || (id > 3) ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    else {
        return ptrfunc[id](nrhs, spm, x0, ldx0, b, ldb, x, ldx );
    }
}

/**
 *******************************************************************************
 *
 * @brief Scale the spm.
 *
 * A = alpha * A
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *           The scaling parameter.
 *
 * @param[inout] spm
 *          The sparse matrix to scal.
 *
 *******************************************************************************/
void
spmScalMatrix(const double alpha, pastix_spm_t* spm)
{
    switch(spm->flttype)
    {
    case PastixPattern:
        break;
    case PastixFloat:
        s_spmScal((float)alpha, spm);
        break;
    case PastixComplex32:
        c_spmScal((float)alpha, spm);
        break;
    case PastixComplex64:
        z_spmScal(alpha, spm);
        break;
    case PastixDouble:
    default:
        d_spmScal(alpha, spm);
    }
}

/**
 *******************************************************************************
 *
 * @brief Scale a vector according to the spm type.
 *
 * x = alpha * x
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *           The scaling parameter.
 *
 * @param[in] spm
 *          The spm structure to know the type of the vector.
 *
 * @param[inout] x
 *          The vector to scal.
 *
 *******************************************************************************/
void
spmScalVector(const double alpha, pastix_spm_t* spm, void *x)
{
    switch(spm->flttype)
    {
    case PastixPattern:
        break;
    case PastixFloat:
        cblas_sscal(spm->n, (float)alpha, x, 1);
        break;
    case PastixComplex32:
        cblas_csscal(spm->n, (float)alpha, x, 1);
        break;
    case PastixComplex64:
        cblas_zdscal(spm->n, alpha, x, 1);
        break;
    case PastixDouble:
    default:
        cblas_dscal(spm->n, alpha, x, 1);
    }
}

/**
 * @}
 */
