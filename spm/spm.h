/**
 *
 * @file spm.h
 *
 *  PaStiX sparse matrix routines to handle different format of sparse matrices.
 *  $COPYRIGHTS$
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _SPM_H_
#define _SPM_H_

/**
 * @ingroup pastix_spm
 *
 * @brief The list of matrix driver reader and generators
 *
 */
typedef enum pastix_driver_e {
    PastixDriverRSA, /* ok */
    PastixDriverCCC,//
    PastixDriverRCC,//
    PastixDriverOlaf,//
    PastixDriverPeer,//
    PastixDriverHB, /* ok */
    PastixDriverIJV, /* ok */
    PastixDriverMM, /* ok */
    PastixDriverDMM, /* ok */
    PastixDriverPetscS, /* ok */
    PastixDriverPetscU, /* ok */
    PastixDriverPetscH, /* ok */
    PastixDriverCSCD,//
    PastixDriverLaplacian, /* ok */
    PastixDriverXLaplacian, /* ok */
    PastixDriverBRGM,//
    PastixDriverBRGMD,//
    PastixDriverGraph
} pastix_driver_t;

/**
 * @ingroup pastix_spm
 *
 * @brief The sparse matrix data structure
 *
 * This structure describes matrices with different characteristics that can be useful to any solver:
 *     - the storage format (PastixCSC, PastixCSR or PastixIJV)
 *     - the properties (PastixGeneral, PastixHermitian, PastixSymmetric)
 *     - the base value (0 in C, 1 in Fortran for examples)
 *
 * It is also possible to describe a matrix with constant or varaibel degrees of freedom.
 *
 */
struct pastix_spm_s {
    int               mtxtype; /**< Matrix structure: PastixGeneral, PastixSymmetric
                                    or PastixHermitian.                                            */
    pastix_coeftype_t flttype; /**< avals datatype: PastixPattern, PastixFloat, PastixDouble,
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
    pastix_order_t    layout;  /**< PastixColMajor, or PastixRowMajor                              */

    pastix_int_t     *colptr;  /**< List of indirections to rows for each vertex [+baseval]        */
    pastix_int_t     *rowptr;  /**< List of edges for each vertex [+baseval]                       */
    pastix_int_t     *loc2glob;/**< Corresponding numbering from local to global [+baseval]        */
    void             *values;  /**< Values stored in the matrix                                    */
};

/**
 * Integer arrays subroutines
 */
pastix_int_t *spmIntConvert( pastix_int_t n, int *input );

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Sorts in ascending order array of element composed of one single
 * pastix_int_t with a single key value.
 *
 *******************************************************************************
 *
 * @param[in,out] pbase TODO
 *
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 *******************************************************************************
 */
void          spmIntSort1Asc1(void * const pbase, const pastix_int_t n);

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Sorts in ascending order array of element composed of two
 * pastix_int_t by ascending order. The first value is used as key.
 *
 *******************************************************************************
 *
 * @param[in,out] pbase TODO
 *
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 *******************************************************************************
 */
void          spmIntSort2Asc1(void * const pbase, const pastix_int_t n);

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Sorts in ascending order array of element composed of two
 * pastix_int_t by ascending order. Both values are used as key.
 *
 *******************************************************************************
 *
 * @param[in,out] pbase TODO
 *
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 *******************************************************************************
 */
void          spmIntSort2Asc2(void * const pbase, const pastix_int_t n);

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_dev
 *
 * @brief Subroutines to print elements of spm structures
 *
 *******************************************************************************
 * @param[in] f Pointer to the file
 *
 * @param[in] i TODO
 *
 * @param[in] j TODO
 *
 * @param[in] A TODO
 *******************************************************************************
 */
static inline void z_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, pastix_complex64_t A ){
    fprintf( f, "%ld %ld %e %e\n", (long)i, (long)j, creal(A), cimag(A) );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_dev
 *
 * @brief Subroutines to print elements of spm structures
 *
 *******************************************************************************
 * @param[in] f Pointer to the file
 *
 * @param[in] i TODO
 *
 * @param[in] j TODO
 *
 * @param[in] A TODO
 *******************************************************************************
 */
static inline void c_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, pastix_complex32_t A ){
    fprintf( f, "%ld %ld %e %e\n", (long)i, (long)j, crealf(A), cimagf(A) );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_dev
 *
 * @brief Subroutines to print elements of spm structures
 *
 *******************************************************************************
 * @param[in] f Pointer to the file
 *
 * @param[in] i TODO
 *
 * @param[in] j TODO
 *
 * @param[in] A TODO
 *******************************************************************************
 */
static inline void d_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, double A ){
    fprintf( f, "%ld %ld %e\n", (long)i, (long)j, A );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_dev
 *
 * @brief Subroutines to print elements of spm structures
 *
 *******************************************************************************
 * @param[in] f Pointer to the file
 *
 * @param[in] i TODO
 *
 * @param[in] j TODO
 *
 * @param[in] A TODO
 *******************************************************************************
 */
static inline void s_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, float A ){
    fprintf( f, "%ld %ld %e\n", (long)i, (long)j, A );
}

/**
 * SPM subroutines
 */
int spmLoad( pastix_spm_t *spm, FILE *infile );
int spmSave( pastix_spm_t *spm, FILE *outfile );

void spmPrint( FILE *f, const pastix_spm_t *spm );

int spmGenRHS(int type, int nrhs, const pastix_spm_t *spm, void *x, int ldx, void *b, int ldb );
int spmCheckAxb( int nrhs, const pastix_spm_t *spm, void *x0, int ldx0, void *b, int ldb, const void *x, int ldx );

void          spmInit( pastix_spm_t *spm );
void          spmExit( pastix_spm_t *spm );
pastix_spm_t *spmCopy( const pastix_spm_t *spm );
void          spmBase( pastix_spm_t *spm, int baseval );
void          spmUpdateComputedFields( pastix_spm_t *spm );
int           spmConvert( int ofmttype, pastix_spm_t *ospm );
void *        spm2Dense( const pastix_spm_t *spm );
pastix_int_t  spmFindBase( const pastix_spm_t *spm );
double        spmNorm( int ntype, const pastix_spm_t *spm );
int           spmMatVec(const pastix_trans_t trans, const void *alpha, const pastix_spm_t *spm, const void *x, const void *beta, void *y );
void          spmScal( const pastix_complex64_t alpha, pastix_spm_t* spm );

int           spmSort( pastix_spm_t *spm );
pastix_int_t  spmMergeDuplicate( pastix_spm_t *spm );
pastix_int_t  spmSymmetrize( pastix_spm_t *spm );

pastix_spm_t *spmCheckAndCorrect( pastix_spm_t *spm );

/**
 * SPM multi-dof subroutines
 */
pastix_spm_t *spmExpand(const pastix_spm_t* spm);
pastix_spm_t *spmDofExtend( const int type, const int dof, const pastix_spm_t *spm );

/**
 * SPM read driver access
 */
int spmReadDriver( pastix_driver_t  driver,
                   char            *filename,
                   pastix_spm_t    *spm,
                   MPI_Comm         pastix_comm );

#endif /* _SPM_H_ */
