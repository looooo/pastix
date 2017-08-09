/*
 *  File: pastix/datatypes.h
 *
 *  Definitions of the datatypes used in PaStiX
 *
 *  Authors:
 *    Mathieu  Faverge    - faverge@labri.fr
 *    Xavier   LACOSTE    - lacoste@labri.fr
 *    Pierre   RAMET      - ramet@labri.fr
 *
 *  Dates:
 *    Version 0.0 - from 08 may 1998
 *                  to   08 jan 2001
 *    Version 1.0 - from 06 jun 2002
 *                  to   06 jun 2002
 */

#ifndef _PASTIX_DATATYPES_H_
#define _PASTIX_DATATYPES_H_
#include <inttypes.h>

/** ****************************************************************************
 * Integers
 **/
#if defined(PASTIX_INT64)

typedef int64_t  pastix_int_t;
typedef uint64_t pastix_uint_t;
#define PASTIX_MPI_INT MPI_INTEGER8

#elif defined(PASTIX_INT32)

typedef int32_t  pastix_int_t;
typedef uint32_t pastix_uint_t;
#define PASTIX_MPI_INT MPI_INTEGER4

#elif defined(PASTIX_LONG)

typedef long          pastix_int_t;
typedef unsigned long pastix_uint_t;
#define PASTIX_MPI_INT MPI_LONG

#else

typedef int          pastix_int_t;
typedef unsigned int pastix_uint_t;
#define PASTIX_MPI_INT MPI_INT

#endif

#if !defined(INTSIZEBITS)
#  define INTSIZEBITS   (sizeof (pastix_int_t) << 3)
#endif

#if !defined(INTVALMAX)
#  define INTVALMAX     ((pastix_int_t) (((pastix_uint_t) 1 << (INTSIZEBITS - 1)) - 1))
#endif

/** ****************************************************************************
 * Double that are not converted through precision generator functions
 **/
typedef double pastix_fixdbl_t;

/** ****************************************************************************
 * Complex numbers (Extracted from PaRSEC project)
 **/
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
/* Windows and non-Intel compiler */
#include <complex>
typedef std::complex<float>  pastix_complex32_t;
typedef std::complex<double> pastix_complex64_t;
#else
typedef float  _Complex      pastix_complex32_t;
typedef double _Complex      pastix_complex64_t;
#endif

#if !defined(__cplusplus) && defined(HAVE_COMPLEX_H)
#include <complex.h>
#else

#ifdef __cplusplus
extern "C" {
#endif

/* These declarations will not clash with what C++ provides because
 * the names in C++ are name-mangled. */

extern double cabs     (pastix_complex64_t z);
extern double creal    (pastix_complex64_t z);
extern double cimag    (pastix_complex64_t z);

extern float  cabsf    (pastix_complex32_t z);
extern float  crealf   (pastix_complex32_t z);
extern float  cimagf   (pastix_complex32_t z);

extern pastix_complex64_t conj  (pastix_complex64_t z);
extern pastix_complex64_t csqrt (pastix_complex64_t z);

extern pastix_complex32_t conjf (pastix_complex32_t z);
extern pastix_complex32_t csqrtf(pastix_complex32_t z);

#ifdef __cplusplus
}
#endif

#endif /* HAVE_COMPLEX_H */


static inline size_t
pastix_size_of(pastix_coeftype_t type)
{
    switch(type) {
    case PastixFloat:     return   sizeof(float);
    case PastixDouble:    return   sizeof(double);
    case PastixComplex32: return 2*sizeof(float);
    case PastixComplex64: return 2*sizeof(double);
    default:
        fprintf(stderr, "pastix_size_of: invalid type parameter\n");
        assert(0);
        return 0;
    }
}

/** ****************************************************************************
 * Pastix data structures
 **/

/* Sparse matrix */
struct pastix_spm_s;
typedef struct pastix_spm_s pastix_spm_t;

/* Main structure of the pastix solver associated to a given problem */
struct pastix_data_s;
typedef struct pastix_data_s pastix_data_t;

/* Graph structure (No values) */
struct pastix_graph_s;
typedef struct pastix_graph_s pastix_graph_t;

/* Ordering structure */
struct pastix_order_s;
typedef struct pastix_order_s pastix_order_t;

/* Solver matrix structure to store L(U)*/
struct solver_matrix_s;
typedef struct solver_matrix_s SolverMatrix;

#endif /* _PASTIX_DATATYPES_H_ */
