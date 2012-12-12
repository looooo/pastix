/**
 *
 * @file plasma.h
 *
 *  PLASMA main header
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_H_
#define _PLASMA_H_

#define PLASMA_VERSION_MAJOR 2
#define PLASMA_VERSION_MINOR 4
#define PLASMA_VERSION_MICRO 6

/* CBLAS requires for scalar arguments to be passed by address rather than by value */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( _val_ ) &(_val_)
#endif

/** ****************************************************************************
 *  PLASMA constants - precisions
 **/
#define PlasmaByte          0
#define PlasmaInteger       1
#define PlasmaRealFloat     2
#define PlasmaRealDouble    3
#define PlasmaComplexFloat  4
#define PlasmaComplexDouble 5

/** ****************************************************************************
 *  PLASMA types
 **/
typedef int  PLASMA_enum;
typedef int  PLASMA_bool;
typedef long PLASMA_index;
typedef long PLASMA_size;

/** ****************************************************************************
 * PLASMA Complex numbers
 **/
#define PLASMA_HAS_COMPLEX_H 1

#if defined(_WIN32) 
# include <float.h>
# if defined(__INTEL_COMPILER) 
    /* Fix name conflict within the cabs prototype (_Complex) that
     * conflicts with a C99 keyword.  */
    #define _Complex __ConflictingComplex
    #include <math.h>
    #undef _Complex 
    #undef complex
    typedef float  _Complex PLASMA_Complex32_t;
    typedef double _Complex PLASMA_Complex64_t;
# else 
    /* Use MS VC complex class */
    #include <complex>
    typedef std::complex<float> PLASMA_Complex32_t;
    typedef std::complex<double> PLASMA_Complex64_t;
    /* For LAPACKE lapacke.h force usage of Windows C++ Complex types */
    #define LAPACK_COMPLEX_CUSTOM
    #define lapack_complex_float std::complex<float>
    #define lapack_complex_double std::complex<double>
    #undef PLASMA_HAS_COMPLEX_H 
# endif
# define isnan _isnan
# define isinf !_finite

#else /* _WIN32 */

    typedef float  _Complex PLASMA_Complex32_t;
    typedef double _Complex PLASMA_Complex64_t;

#endif 

/* Sun doesn't ship the complex.h header. Sun Studio doesn't have it and older GCC compilers don't have it either. */
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC) || defined(sun) || defined(__sun)
#undef PLASMA_HAS_COMPLEX_H
#endif

#if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1)) 
#define PLASMA_DEPRECATED  __attribute__((__deprecated__))
#else
#define PLASMA_DEPRECATED
#endif /* __GNUC__ */

#ifdef PLASMA_HAS_COMPLEX_H
#include <complex.h>
#else

#ifdef __cplusplus
extern "C" {
#endif

/* These declarations will not clash with what C++ provides because the names in C++ are name-mangled. */
#if !defined(_WIN32) 
extern double cabs(PLASMA_Complex64_t z);
extern PLASMA_Complex64_t conj(PLASMA_Complex64_t z);
#endif
extern float cabsf(PLASMA_Complex32_t z);
extern double cimag(PLASMA_Complex64_t z);
extern double creal(PLASMA_Complex64_t z);

#ifdef __cplusplus
}
#endif

#endif

#include <quark.h>

/** ****************************************************************************
 *  Tile matrix descriptor
 *
 *  Matrices are stored in a contiguous data chunk containning in order
 *  A11, A21, A12, A22 with :
 *
 *           n1      n2
 *      +----------+---+
 *      |          |   |    With m1 = lm - (lm%mb)
 *      |          |   |         m2 = lm%mb
 *  m1  |    A11   |A12|         n1 = ln - (ln%mb)
 *      |          |   |         n2 = ln%nb
 *      |          |   |
 *      +----------+---+
 *  m2  |    A21   |A22|
 *      +----------+---+
 *
 */
typedef struct plasma_desc_t {
    void *mat;          // pointer to the beginning of the matrix
    size_t A21;        // pointer to the beginning of the matrix A21
    size_t A12;        // pointer to the beginning of the matrix A12
    size_t A22;        // pointer to the beginning of the matrix A22
    PLASMA_enum dtyp;   // precision of the matrix
    int mb;             // number of rows in a tile
    int nb;             // number of columns in a tile
    int bsiz;           // size in elements including padding
    int lm;             // number of rows of the entire matrix
    int ln;             // number of columns of the entire matrix
    int lm1;            // number of tile rows of the A11 matrix - derived parameter
    int ln1;            // number of tile columns of the A11 matrix - derived parameter
    int lmt;            // number of tile rows of the entire matrix - derived parameter
    int lnt;            // number of tile columns of the entire matrix - derived parameter
    int i;              // row index to the beginning of the submatrix
    int j;              // column index to the beginning of the submatrix
    int m;              // number of rows of the submatrix
    int n;              // number of columns of the submatrix
    int mt;             // number of tile rows of the submatrix - derived parameter
    int nt;             // number of tile columns of the submatrix - derived parameter
} PLASMA_desc;

/** ****************************************************************************
 *  PLASMA request uniquely identifies each asynchronous function call.
 **/
typedef struct plasma_request_t {
    PLASMA_enum status; // PLASMA_SUCCESS or appropriate error code
} PLASMA_request;

#define PLASMA_REQUEST_INITIALIZER {PLASMA_SUCCESS}

/** ****************************************************************************
 *  PLASMA sequence uniquely identifies a set of asynchronous function calls
 *  sharing common exception handling.
 **/
typedef struct plasma_sequence_t {
    Quark_Sequence *quark_sequence; // QUARK sequence associated with PLASMA sequence
    PLASMA_bool     status;         // PLASMA_SUCCESS or appropriate error code
    PLASMA_request *request;        // failed request
} PLASMA_sequence;

/** ****************************************************************************
 *
 *  PLASMA constants - CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz),
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/).
 *
 **/
#define PlasmaRM            101
#define PlasmaCM            102
#define PlasmaCCRB          103
#define PlasmaCRRB          104
#define PlasmaRCRB          105
#define PlasmaRRRB          106

#define PlasmaNoTrans       111
#define PlasmaTrans         112
#define PlasmaConjTrans     113

#define PlasmaUpper         121
#define PlasmaLower         122
#define PlasmaUpperLower    123

#define PlasmaNonUnit       131
#define PlasmaUnit          132

#define PlasmaLeft          141
#define PlasmaRight         142

#define PlasmaOneNorm       171
#define PlasmaRealOneNorm   172
#define PlasmaTwoNorm       173
#define PlasmaFrobeniusNorm 174
#define PlasmaInfNorm       175
#define PlasmaRealInfNorm   176
#define PlasmaMaxNorm       177
#define PlasmaRealMaxNorm   178

#define PlasmaDistUniform   201
#define PlasmaDistSymmetric 202
#define PlasmaDistNormal    203

#define PlasmaHermGeev      241
#define PlasmaHermPoev      242
#define PlasmaNonsymPosv    243
#define PlasmaSymPosv       244

#define PlasmaNoPacking     291
#define PlasmaPackSubdiag   292
#define PlasmaPackSupdiag   293
#define PlasmaPackColumn    294
#define PlasmaPackRow       295
#define PlasmaPackLowerBand 296
#define PlasmaPackUpeprBand 297
#define PlasmaPackAll       298

#define PlasmaNoVec         301
#define PlasmaVec           302
#define PlasmaIvec          303

#define PlasmaForward       391
#define PlasmaBackward      392

#define PlasmaColumnwise    401
#define PlasmaRowwise       402
#define PlasmaTrd          1001
#define PlasmaBrd          1002

#define PlasmaW             501
#define PlasmaA2            502

#define plasma_const_neg(const) (((const-1)^0x01)+1)

/** ****************************************************************************
 *  PLASMA constants - boolean
 **/
#define PLASMA_FALSE 0
#define PLASMA_TRUE  1

/** ****************************************************************************
 *  State machine switches
 **/
#define PLASMA_WARNINGS   1
#define PLASMA_ERRORS     2
#define PLASMA_AUTOTUNING 3
#define PLASMA_DAG        4

/** ****************************************************************************
 *  PLASMA constants - configuration parameters
 **/
#define PLASMA_CONCURRENCY      1
#define PLASMA_TILE_SIZE        2
#define PLASMA_INNER_BLOCK_SIZE 3
#define PLASMA_SCHEDULING_MODE  4
#define PLASMA_HOUSEHOLDER_MODE 5
#define PLASMA_HOUSEHOLDER_SIZE 6
#define PLASMA_TRANSLATION_MODE 7

#define PLASMA_STATIC_SCHEDULING  1
#define PLASMA_DYNAMIC_SCHEDULING 2

#define PLASMA_FLAT_HOUSEHOLDER 1
#define PLASMA_TREE_HOUSEHOLDER 2

#define PLASMA_INPLACE    1
#define PLASMA_OUTOFPLACE 2

/** ****************************************************************************
 *  PLASMA constants - success & error codes
 **/
#define PLASMA_SUCCESS                 0
#define PLASMA_ERR_NOT_INITIALIZED  -101
#define PLASMA_ERR_REINITIALIZED    -102
#define PLASMA_ERR_NOT_SUPPORTED    -103
#define PLASMA_ERR_ILLEGAL_VALUE    -104
#define PLASMA_ERR_NOT_FOUND        -105
#define PLASMA_ERR_OUT_OF_RESOURCES -106
#define PLASMA_ERR_INTERNAL_LIMIT   -107
#define PLASMA_ERR_UNALLOCATED      -108
#define PLASMA_ERR_FILESYSTEM       -109
#define PLASMA_ERR_UNEXPECTED       -110
#define PLASMA_ERR_SEQUENCE_FLUSHED -111

/** ****************************************************************************
 *  Math function prototypes
 **/
#include <plasma_z.h>
#include <plasma_d.h>
#include <plasma_c.h>
#include <plasma_s.h>
#include <plasma_zc.h>
#include <plasma_ds.h>

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Auxiliary function prototypes
 **/
int PLASMA_Version(int *ver_major, int *ver_minor, int *ver_micro);
int PLASMA_Enable(PLASMA_enum lever);
int PLASMA_Disable(PLASMA_enum lever);
int PLASMA_Set(PLASMA_enum param, int value);
int PLASMA_Get(PLASMA_enum param, int *value);
int PLASMA_Init(int cores);
int PLASMA_Init_Affinity(int cores, int *bindtab);
int PLASMA_Finalize();
int PLASMA_Desc_Create(PLASMA_desc **desc, void *mat, PLASMA_enum dtyp, int mb, int nb, int bsiz, int lm, int ln, int i, int j, int m, int n);
int PLASMA_Desc_Destroy(PLASMA_desc **desc);
int PLASMA_Lapack_to_Tile(void *Af77, int LDA, PLASMA_desc *A);
int PLASMA_Tile_to_Lapack(PLASMA_desc *A, void *Af77, int LDA);

/** ****************************************************************************
 *  Workspace deallocation
 **/
int PLASMA_Dealloc_Handle(void **handle);
int PLASMA_Dealloc_Handle_Tile(PLASMA_desc **desc);

/** ****************************************************************************
 *  Sequence
 **/
int PLASMA_Sequence_Create(PLASMA_sequence **sequence);
int PLASMA_Sequence_Destroy(PLASMA_sequence *sequence);
int PLASMA_Sequence_Wait(PLASMA_sequence *sequence);
int PLASMA_Sequence_Flush(PLASMA_sequence *sequence, PLASMA_request *request);

#ifdef __cplusplus
}
#endif

#endif
