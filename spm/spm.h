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
 * @enum pastix_driver_e
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
 * @struct pastix_spm_s
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
    int                mtxtype;/**< Matrix structure: PastixGeneral, PastixSymmetric
                                    or PastixHermitian.                                         */
    pastix_coeftype_t flttype; /**< avals datatype: PastixPattern, PastixFloat, PastixDouble,
                                    PastixComplex32 or PastixComplex64                          */
    pastix_fmttype_t  fmttype; /**< Matrix storage format: PastixCSC, PastixCSR, PastixIJV      */

    pastix_int_t      gN;      /**< Global number of vertices in the compressed graph           */
    pastix_int_t      n;       /**< Local number of vertices in the compressed graph            */
    pastix_int_t      gnnz;    /**< Global number of non zeroes in the compressed graph         */
    pastix_int_t      nnz;     /**< Local number of non zeroes in the compressed graph          */

    pastix_int_t      gNexp;   /**< Global number of vertices in the compressed graph           */
    pastix_int_t      nexp;    /**< Local number of vertices in the compressed graph            */
    pastix_int_t      gnnzexp; /**< Global number of non zeroes in the compressed graph         */
    pastix_int_t      nnzexp;  /**< Local number of non zeroes in the compressed graph          */

    pastix_int_t      dof;     /**< Number of degrees of freedom per unknown,
                                    if > 0, constant degree of freedom
                                    otherwise, irregular degree of freedom (refer to dofs)      */
    pastix_int_t     *dofs;    /**< Number of degrees of freedom per unknown (NULL, if dof > 0) */
    pastix_order_t    layout;  /**< PastixColMajor, or PastixRowMajor                           */

    pastix_int_t     *colptr;  /**< List of indirections to rows for each vertex                */
    pastix_int_t     *rowptr;  /**< List of edges for each vertex                               */
    pastix_int_t     *loc2glob;/**< Corresponding numbering from local to global                */
    void             *values;  /**< Values stored in the matrix                                 */
};

int
csc_load( pastix_int_t  *n,
          pastix_int_t **colptr,
          pastix_int_t **rows,
          int           *valtype,
          void         **values,
          int           *dof,
          FILE          *infile );

int
csc_save( pastix_int_t  n,
          pastix_int_t *colptr,
          pastix_int_t *rows,
          int           ft,
          void         *values,
          int           dof,
          FILE         *outfile );
/**
 * Integer arrays subroutines
 */
pastix_int_t *spmIntConvert( pastix_int_t n, int *input );
void          spmIntSort1Asc1(void * const pbase, const pastix_int_t n);
void          spmIntSort2Asc1(void * const pbase, const pastix_int_t n);
void          spmIntSort2Asc2(void * const pbase, const pastix_int_t n);


int spmLoad( pastix_spm_t *spm, FILE *infile );
int spmSave( pastix_spm_t *spm, FILE *outfile );

int spmGenRHS(int type, int nrhs, const pastix_spm_t *spm, void *x, int ldx, void *b, int ldb );
int spmCheckAxb( int nrhs, const pastix_spm_t *spm, void *x0, int ldx0, void *b, int ldb, const void *x, int ldx );

void          spmInit( pastix_spm_t *spm );
void          spmExit( pastix_spm_t *spm );
pastix_spm_t *spmCopy( const pastix_spm_t *spm );
void          spmBase( pastix_spm_t *spm, int baseval );
int           spmConvert( int ofmttype, pastix_spm_t *ospm );
pastix_int_t  spmFindBase( const pastix_spm_t *spm );
double        spmNorm( int ntype, const pastix_spm_t *spm );
int           spmMatVec(const pastix_trans_t trans, const void *alpha, const pastix_spm_t *spm, const void *x, const void *beta, void *y );

int           spmSort( pastix_spm_t *spm );
pastix_int_t  spmMergeDuplicate( pastix_spm_t *spm );
pastix_int_t  spmSymmetrize( pastix_spm_t *spm );

pastix_spm_t *spmCheckAndCorrect( pastix_spm_t *spm );

pastix_spm_t *spmExpand(const pastix_spm_t* spm);

int spmReadDriver( pastix_driver_t  driver,
                   char            *filename,
                   pastix_spm_t    *spm,
                   MPI_Comm         pastix_comm );

#endif /* _SPM_H_ */
