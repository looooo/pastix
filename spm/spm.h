/**
 *
 * @file spm.h
 *
 * SParse Matrix package header.
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
 *   @brief Describe all the internals routines of the SParse Matrix package.
 *
 *   This library provides a set of subroutines to manipulate sparse matrices in
 *   different format such as compressed sparse column (CSC), compressed sparse
 *   row (CSR), or coordinate (IJV) with single or multiple degrees of freedom
 *   per unknown. It provides basic BLAS 1 and BLAS 2 functions for those
 *   matrices, as well as norms computations and converter tools.
 *
 **/
#ifndef _spm_h_
#define _spm_h_

#include "pastix/api.h"

/**
 *
 * @brief The sparse matrix data structure
 *
 * This structure describes matrices with different characteristics that can be useful to any solver:
 *     - the storage format (PastixCSC, PastixCSR or PastixIJV)
 *     - the properties (PastixGeneral, PastixHermitian, PastixSymmetric)
 *     - the base value (0 in C or 1 in Fortran)
 *
 * It is also possible to describe a matrix with constant or variable degrees of freedom.
 *
 */
typedef struct pastix_spm_s {
    pastix_mtxtype_t  mtxtype; /**< Matrix structure: PastixGeneral, PastixSymmetric
                                    or PastixHermitian.                                            */
    pastix_coeftype_t flttype; /**< values datatype: PastixPattern, PastixFloat, PastixDouble,
                                    PastixComplex32 or PastixComplex64                             */
    pastix_fmttype_t  fmttype; /**< Matrix storage format: PastixCSC, PastixCSR, PastixIJV         */

    pastix_int_t      gN;      /**< Global number of vertices in the compressed graph (Computed)   */
    pastix_int_t      n;       /**< Local number of vertices in the compressed graph               */
    pastix_int_t      gnnz;    /**< Global number of non zeroes in the compressed graph (Computed) */
    pastix_int_t      nnz;     /**< Local number of non zeroes in the compressed graph             */

    pastix_int_t      gNexp;   /**< Global number of vertices in the compressed graph (Computed)   */
    pastix_int_t      nexp;    /**< Local number of vertices in the compressed graph (Computed)    */
    pastix_int_t      gnnzexp; /**< Global number of non zeroes in the compressed graph (Computed) */
    pastix_int_t      nnzexp;  /**< Local number of non zeroes in the compressed graph (Computed)  */

    pastix_int_t      dof;     /**< Number of degrees of freedom per unknown,
                                    if > 0, constant degree of freedom
                                    otherwise, irregular degree of freedom (refer to dofs)         */
    pastix_int_t     *dofs;    /**< Array of the first column of each element in the
                                    expanded matrix [+baseval]                                     */
    pastix_layout_t   layout;  /**< PastixColMajor, or PastixRowMajor                              */

    pastix_int_t     *colptr;  /**< List of indirections to rows for each vertex [+baseval]        */
    pastix_int_t     *rowptr;  /**< List of edges for each vertex [+baseval]                       */
    pastix_int_t     *loc2glob;/**< Corresponding numbering from local to global [+baseval]        */
    void             *values;  /**< Values stored in the matrix                                    */
} pastix_spm_t;

/**
 * @name SPM basic subroutines
 * @{
 */
void          spmInit( pastix_spm_t *spm );
void          spmExit( pastix_spm_t *spm );

pastix_spm_t *spmCopy( const pastix_spm_t *spm );
void          spmBase( pastix_spm_t *spm, int baseval );
pastix_int_t  spmFindBase( const pastix_spm_t *spm );
int           spmConvert( int ofmttype, pastix_spm_t *ospm );
void          spmUpdateComputedFields( pastix_spm_t *spm );
void          spmGenFakeValues( pastix_spm_t *spm );

/**
 * @}
 * @name SPM BLAS subroutines
 * @{
 */
double        spmNorm( pastix_normtype_t ntype, const pastix_spm_t *spm );
int           spmMatVec( pastix_trans_t trans, const void *alpha, const pastix_spm_t *spm, const void *x, const void *beta, void *y );
void          spmScalMatrix( double alpha, pastix_spm_t *spm );
void          spmScalVector( double alpha, pastix_spm_t *spm, void *x );
void          spmScalRHS( pastix_coeftype_t flt, double alpha, pastix_int_t m, pastix_int_t n, void *A, pastix_int_t lda );

/**
 * @}
 * @name SPM subroutines to check format
 * @{
 */
int           spmSort( pastix_spm_t *spm );
pastix_int_t  spmMergeDuplicate( pastix_spm_t *spm );
pastix_int_t  spmSymmetrize( pastix_spm_t *spm );
pastix_spm_t *spmCheckAndCorrect( pastix_spm_t *spm );

/**
 * @}
 * @name SPM subroutines to check factorization/solve
 * @{
 */
int           spmGenRHS( pastix_rhstype_t type, pastix_int_t nrhs, const pastix_spm_t *spm, void *x, pastix_int_t ldx, void *b, pastix_int_t ldb );
int           spmCheckAxb( pastix_int_t nrhs, const pastix_spm_t *spm, void *x0, pastix_int_t ldx0, void *b, pastix_int_t ldb, const void *x, pastix_int_t ldx );

/**
 * @}
 * @name SPM subroutines to manipulate integers arrays
 * @{
 */
pastix_int_t *spmIntConvert(   pastix_int_t n, int *input );
void          spmIntSort1Asc1( void * const pbase, const pastix_int_t n );
void          spmIntSort2Asc1( void * const pbase, const pastix_int_t n );
void          spmIntSort2Asc2( void * const pbase, const pastix_int_t n );

/**
 * @}
 * @name SPM IO subroutines
 * @{
 */
int           spmLoad(       pastix_spm_t *spm, FILE *infile );
int           spmSave( const pastix_spm_t *spm, FILE *outfile );

/**
 * @}
 * @name SPM driver
 * @{
 */
int           spmReadDriver( pastix_driver_t  driver,
                             const char      *filename,
                             pastix_spm_t    *spm,
                             MPI_Comm         pastix_comm );
/**
 * @}
 * @name SPM debug subroutines
 * @{
 */
void *        spm2Dense   ( const pastix_spm_t *spm );
void          spmPrint    ( const pastix_spm_t *spm, FILE *f );
void          spmPrintInfo( const pastix_spm_t *spm, FILE *f );
pastix_spm_t *spmExpand   ( const pastix_spm_t *spm );
pastix_spm_t *spmDofExtend( const pastix_spm_t *spm, const int type, const int dof );

/**
 * @}
 */

/**
 * @}
 */

/**
 * @name SPM dev printing subroutines
 * @{
 *
 */

/**
 * @ingroup spm_dev_print
 * @brief Subroutines to print one element of an spm structure
 *
 * @param[in] f Pointer to the file
 * @param[in] i Row index of the element
 * @param[in] j Column index of the element
 * @param[in] A Value of the element A|i,j]
 *
 * Double complex case
 *
 */
static inline void z_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, pastix_complex64_t A ){
    fprintf( f, "%ld %ld %e %e\n", (long)i, (long)j, creal(A), cimag(A) );
}

/**
 * @copydoc z_spmPrintElt
 * @details Single complex case
 */
static inline void c_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, pastix_complex32_t A ){
    fprintf( f, "%ld %ld %e %e\n", (long)i, (long)j, crealf(A), cimagf(A) );
}
/**
 * @copydoc z_spmPrintElt
 * @details Double real case
 */
static inline void d_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, double A ){
    fprintf( f, "%ld %ld %e\n", (long)i, (long)j, A );
}
/**
 * @copydoc z_spmPrintElt
 * @details Single real case
 */
static inline void s_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, float A ){
    fprintf( f, "%ld %ld %e\n", (long)i, (long)j, A );
}
/**
 * @copydoc z_spmPrintElt
 * @details Pattern case
 *
 * @remark: uses a macro to avoid accessing A that would generate segfault.
 */
#define p_spmPrintElt( f, i, j, A ) {                           \
        fprintf( f, "%ld %ld\n", (long)(i), (long)(j) );        \
    }

/**
 * @}
 */
#endif /* _spm_h_ */
