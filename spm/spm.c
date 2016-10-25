/**
 *
 * @file spm.c
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

#include "z_spm.h"
#include "c_spm.h"
#include "d_spm.h"
#include "s_spm.h"
#include "p_spm.h"

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


/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Init the spm structure given as parameter
 *
 *******************************************************************************
 *
 * @param[in,out] spm
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

    spm->dof       = 1;
    spm->dofs      = NULL;
    spm->colmajor  = 1;

    spm->colptr   = NULL;
    spm->rowptr   = NULL;
    spm->loc2glob = NULL;
    spm->values   = NULL;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Free the spm structure
 *
 *******************************************************************************
 *
 * @param[in,out] spm
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
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Rebase the spm
 *
 * Rebase the arrays of the spm to the given value. If the value is equal to the
 * original base, then nothing is performed.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
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
        pastix_error_print("spmBase: spm pointer is NULL");
        return;
    }
    if ( (spm->colptr == NULL) ||
         (spm->rowptr == NULL) )
    {
        pastix_error_print("spmBase: spm pointer is not correctly initialized");
        return;
    }
    if ( (baseval != 0) &&
         (baseval != 1) )
    {
        pastix_error_print("spmBase: baseval is incorrect, must be 0 or 1");
        return;
    }

    baseadj = baseval - spmFindBase( spm );
    if (baseadj == 0)
	return;

    n   = spm->n;
    nnz = spm->nnz;

    assert( nnz == (spm->colptr[n] - spm->colptr[0]) );

    for (i = 0; i <= n; i++) {
        spm->colptr[i] += baseadj;
    }
    for (i = 0; i < nnz; i++) {
        spm->rowptr[i] += baseadj;
    }

    if (spm->loc2glob != NULL) {
        for (i = 0; i < n; i++) {
            spm->loc2glob[i] += baseadj;
        }
    }
    return;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
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
 * @ingroup pastix_spm
 *
 * @brief  Convert the storage format of the spm to ofmttype.
 *
 * Convert the storage format of the given sparse matrix from any of the
 * following format: PastixCSC, PastixCSR, or PastixIJV to one of these.
 *
 *******************************************************************************
 *
 * @param[in] ofmttype
 *          The output format of the sparse matrix. It might be PastixCSC,
 *          PastixCSR, or PastixIJV.
 *
 * @param[in,out] spm
 *          The sparse matrix structure to convert.
 *
 ********************************************************************************
 *
 * @return
 *        \retval PASTIX_SUCCESS if the conversion happened successfuly
 *        \retval PASTIX_ERR_BADPARAMETER if one the parameter is incorrect.
 *
 *******************************************************************************/
int
spmConvert( int ofmttype, pastix_spm_t *spm )
{
    if ( spm->dof != 1 ) {
        spm = spmExpand( spm );
    }
    if ( conversionTable[spm->fmttype][ofmttype][spm->flttype] ) {
        return conversionTable[spm->fmttype][ofmttype][spm->flttype]( spm );
    }
    else {
        return PASTIX_SUCCESS;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
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
 *          = PastixMaxNorm: Max norm
 *          = PastixOneNorm: One norm
 *          = PastixInfNorm: Infinity norm
 *          = PastixFrobeniusNorm: Frobenius norm
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 ********************************************************************************
 *
 * @return
 *          \retval the norm described above. Note that for simplicity, even if
 *          the norm of single real or single complex matrix is computed with
 *          single precision, the returned norm is stored in double precision
 *          number.
 *          \retval -1., if the floating point of the sparse matrix is
 *          undefined.
 *
 *******************************************************************************/
double
spmNorm( int ntype,
         const pastix_spm_t *spm )
{
    double tmp;

    if ( spm->dof != 1 ) {
        fprintf(stderr, "WARNING: spm expanded due to non implemented norm for non-expanded spm\n");
        spm = spmExpand( spm );
    }
    switch (spm->flttype) {
    case PastixFloat:
        tmp = (double)s_spmNorm( ntype, spm );
        return tmp;

    case PastixDouble:
        return d_spmNorm( ntype, spm );

    case PastixComplex32:
        tmp = (double)c_spmNorm( ntype, spm );
        return tmp;

    case PastixComplex64:
        return z_spmNorm( ntype, spm );

    default:
        return -1.;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Sort the subarray of edges of each vertex in a CSC or CSR spm
 *
 * This routine sorts the subarray of edges of each vertex in a
 * centralized spm stored in CSC or CSR format. Nothing is performed if IJV
 * format is used.
 *
 * WARNING: This function should NOT be called if dof is greater than 1.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with subarrays of edges sorted by
 *          ascending order.
 *
 ********************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS if the sort was called
 *          \retval PASTIX_ERR_BADPARAMETER, if the given spm was incorrect.
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
 * @ingroup pastix_spm
 *
 * @brief Merge mulitple entries in a spm by summing them.
 *
 * This routine merge the multiple entries in a sparse
 * matrix by suming their values together. The sparse matrix needs to be sorted
 * first (see spmSort()).
 *
 * WARNING: Not implemented for CSR and IJV format.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the reduced sparse matrix of multiple entries were present
 *          in it. The multiple values for a same vertex are sum up together.
 *
 ********************************************************************************
 *
 * @return
 *          \retval If >=0, the number of vertices that were merged
 *          \retval PASTIX_ERR_BADPARAMETER, if the given spm was incorrect.
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
 * @ingroup pastix_spm
 *
 * @brief Symmetrize the pattern of the spm
 *
 * Symmetrize the pattern of the input spm by edges when A(i,j) exists, but
 * A(j,i) does not. When values are associated to the edge, zeroes are added to
 * the values array.
 *
 * WARNING: Not implemented for CSR and IJV format.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the reduced sparse matrix of multiple entries were present
 *          in it. The multiple values for a same vertex are sum up together.
 *
 ********************************************************************************
 *
 * @return
 *          \retval If >=0, the number of vertices that were merged
 *          \retval PASTIX_ERR_BADPARAMETER, if the given spm was incorrect.
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
 * @ingroup pastix_spm
 *
 * @brief Check the correctness of a spm.
 *
 * This routine initializes the sparse matrix to fit the PaStiX requirements. If
 * needed, the format is changed to CSC, the duplicated vertices are merged
 * together by summing their values; the graph is made symmetric for matrices
 * with unsymmetric pattern, new values are set to 0.; Only the lower part is
 * kept for the symmetric matrices.
 *
 * On exit, if no changes have been made, the initial sparse matrix is returned,
 * otherwise a copy of the sparse matrix structured fixed to meet the PaStiX
 * requirements is returned.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The pointer to the sparse matrix structure to check, and correct.
 *          On exit, the subarrays related to each column might have been sorted
 *          by ascending order.
 *
 *******************************************************************************
 *
 * @return
 *          \retval If no modifications were made to the initial matrix
 *                  structure, the one given as parameter is returned
 *          \retval Otherwise, the news sparse matrix structure is returned. It
 *                  must be destroyed with spmExit() and a free of the returned
 *                  pointer.
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
        fprintf(stderr, "spmCheckAndCorrect: %ld entries have been merged\n", (int64_t)count );
    }

    /**
     * If the matrix is symmetric or hermitian, we keep only the upper or lower
     * part, otherwise, we symmetrize the graph to get A+A^t, new values are set
     * to 0.
     */
    if ( newspm->mtxtype == PastixGeneral ) {
        count = spmSymmetrize( newspm );
        if ( count > 0 ) {
            fprintf(stderr, "spmCheckAndCorrect: %ld entries have been added for symmetry\n", (int64_t)count );
        }
    }
    else {
        //spmToLower( newspm );
    }

    /**
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
 * @ingroup pastix_spm
 *
 * @brief Create a copy of the spm
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

    memcpy( newspm, spm, sizeof(pastix_spm_t));

    if(spm->colptr != NULL) {
        newspm->colptr = (pastix_int_t*)malloc((spm->n+1) * sizeof(pastix_int_t));
        memcpy( newspm->colptr, spm->colptr, (spm->n+1) * sizeof(pastix_int_t));
    }
    if(spm->rowptr != NULL) {
        newspm->rowptr = (pastix_int_t*)malloc(spm->nnz * sizeof(pastix_int_t));
        memcpy( newspm->rowptr, spm->rowptr, spm->nnz * sizeof(pastix_int_t));
    }
    if(spm->loc2glob != NULL) {
        newspm->loc2glob = (pastix_int_t*)malloc(spm->n * sizeof(pastix_int_t));
        memcpy( newspm->loc2glob, spm->loc2glob, spm->n * sizeof(pastix_int_t));
    }
    if(spm->values != NULL) {
        size_t valsize = spm->nnzexp * pastix_size_of( spm->flttype );
        newspm->values = malloc(valsize);
        memcpy( newspm->values, spm->values, valsize);
    }
    return newspm;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Expand a multi-dof spm matrix into an spm with constant dof to 1.
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
void
spmExpand(pastix_spm_t* spm)
{
    switch(spm->flttype)
    {
    case PastixFloat:
        s_spmExpand(spm);
        break;
    case PastixComplex32:
        c_spmExpand(spm);
        break;
    case PastixComplex64:
        z_spmExpand(spm);
        break;
    case PastixDouble:
    default:
        d_spmExpand(spm);
        break;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Compute a matrix-vector product
 *
 * Compute the matrix vector product:
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
 *          Specifies whether the matrix spm is transposed, not transposed or
 *          conjugate transposed:
 *          = PastixNoTrans:   A is not transposed;
 *          = PastixTrans:     A is transposed;
 *          = PastixConjTrans: A is conjugate transposed.
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
 * @param[in,out] y
 *          The vector y.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the y vector has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
/**
 * TODO: Maybe we should move down the cast of the parameters to the lowest
 * functions, and simplify this one to have identical calls to all subfunction
 */
int
spmMatVec(      int           trans,
          const void         *alpha,
          const pastix_spm_t *spm,
          const void         *x,
          const void         *beta,
                void         *y )
{
    switch (spm->mtxtype) {
    case PastixHermitian:
        switch (spm->flttype) {
        case PastixFloat:
            return s_spmSyCSCv( *((const float*)alpha), spm, (const float*)x, *((const float*)beta), (float*)y );
        case PastixComplex32:
            return c_spmHeCSCv( *((const pastix_complex32_t*)alpha), spm, (const pastix_complex32_t*)x, *((const pastix_complex32_t*)beta), (pastix_complex32_t*)y );
        case PastixComplex64:
            return z_spmHeCSCv( *((const pastix_complex64_t*)alpha), spm, (const pastix_complex64_t*)x, *((const pastix_complex64_t*)beta), (pastix_complex64_t*)y );
        case PastixDouble:
        default:
            return d_spmSyCSCv( *((const double*)alpha), spm, (const double*)x, *((const double*)beta), (double*)y );
        }
    case PastixSymmetric:
        switch (spm->flttype) {
        case PastixFloat:
            return s_spmSyCSCv( *((const float*)alpha), spm, (const float*)x, *((const float*)beta), (float*)y );
        case PastixComplex32:
            return c_spmSyCSCv( *((const pastix_complex32_t*)alpha), spm, (const pastix_complex32_t*)x, *((const pastix_complex32_t*)beta), (pastix_complex32_t*)y );
        case PastixComplex64:
            return z_spmSyCSCv( *((const pastix_complex64_t*)alpha), spm, (const pastix_complex64_t*)x, *((const pastix_complex64_t*)beta), (pastix_complex64_t*)y );
        case PastixDouble:
        default:
            return d_spmSyCSCv( *((const double*)alpha), spm, (const double*)x, *((const double*)beta), (double*)y );
        }
    case PastixGeneral:
    default:
        switch (spm->flttype) {
        case PastixFloat:
            return s_spmGeCSCv( trans, *((const float*)alpha), spm, (const float*)x, *((const float*)beta), (float*)y );
        case PastixComplex32:
            return c_spmGeCSCv( trans, *((const pastix_complex32_t*)alpha), spm, (const pastix_complex32_t*)x, *((const pastix_complex32_t*)beta), (pastix_complex32_t*)y );
        case PastixComplex64:
            return z_spmGeCSCv( trans, *((const pastix_complex64_t*)alpha), spm, (const pastix_complex64_t*)x, *((const pastix_complex64_t*)beta), (pastix_complex64_t*)y );
        case PastixDouble:
        default:
            return d_spmGeCSCv( trans, *((const double*)alpha), spm, (const double*)x, *((const double*)beta), (double*)y );
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Generate right hand sides.
 *
 * Generate nrhs right hand side vectors associated to a given matrix to test a
 * problem with a solver. The vectors can be initialized randomly, or to get a
 * specific solution.
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
 * @param[in,out] b
 *          b must be an allocated matrix of size at least ldb * nrhs.
 *          On exit, b is initialized as defined by the type parameter.
 *
 * @param[in] ldb
 *          Defines the leading dimension of b when multiple right hand sides
 *          are available. ldb >= spm->n.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the b vector has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmGenRHS( int type, int nrhs,
           const pastix_spm_t  *spm,
           void                *x, int ldx,
           void                *b, int ldb )
{
    static int (*ptrfunc[4])(int, int,
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
 * @ingroup pastix_spm
 *
 * @brief Check backward and forward errors
 *
 * Check the backward error, and the forward error if x0 is provided.
 *
 *******************************************************************************
 *
 * @param[in] nrhs
 *          Defines the number of right hand side that must be generated.
 *
 * @param[in] spm
 *          The sparse matrix uses to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[in,out] x0
 *          If x0 != NULL, the forward error is computed.
 *          On exit, x0 stores (x0-x)
 *
 * @param[in] ldx0
 *          Defines the leading dimension of x0 when multiple right hand sides
 *          are available. ldx0 >= spm->n.
 *
 * @param[in,out] b
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
 * @return
 *      \retval PASTIX_SUCCESS if the b vector has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmCheckAxb( int nrhs,
             const pastix_spm_t  *spm,
                   void *x0, int ldx0,
                   void *b,  int ldb,
             const void *x,  int ldx )
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
