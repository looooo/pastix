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
 * Precisions
 * (Start at 2 for Compatibility with Plasma, just in case)
 **/

#define PastixFloat     2
#define PastixDouble    3
#define PastixComplex32 4
#define PastixComplex64 5

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

#endif /* _PASTIX_DATATYPES_H_ */
