#ifndef PASTIX_RETRO_H
#define PASTIX_RETRO_H

#define pastix_int_t pastix_int_t
#define pastix_float_t double

#define ASSERTDBG( test, mod ) assert(test)
#define ASSERT( test, mod ) assert(test)


/* #define PASTIX_PREFIX_F(x)   _PASTIX_ ## x */
/* #define PASTIX_EXTERN_F(x)   x */
/* #define PASTIX_PREFIX(x)     _PASTIX_ ## x */
/* #define PASTIX_EXTERN(x)     x */

#define CHECK_MPI(call) call
#  define CHECK_THREAD_LEVEL(THREAD_MODE) do {  \
    THREAD_MODE = THREAD_MODE;                  \
  } while (0)

#define CONJ_FLOAT(val) (val)
#define ABS_FLOAT(val) abs(val)

#define BLAS_INT int
#define BLAS_REAL double
#define BLAS_FLOAT double
#define TAB_CHECK_NAN(tab, size)  do {} while (0)

#define COMM_FLOAT MPI_DOUBLE
#define COMM_SUM MPI_SUM
#define PASTIX_MPI_FLOAT MPI_DOUBLE

#define NO_ERR              PASTIX_SUCCESS
#define UNKNOWN_ERR         PASTIX_ERR_UNKNOWN
#define ALLOC_ERR           PASTIX_ERR_ALLOC
#define ASSERT_ERR          PASTIX_ERR_ASSERT
#define NOTIMPLEMENTED_ERR  PASTIX_ERR_NOTIMPLEMENTED
#define OUTOFMEMORY_ERR     PASTIX_ERR_OUTOFMEMORY
#define THREAD_ERR          PASTIX_ERR_THREAD
#define INTERNAL_ERR        PASTIX_ERR_INTERNAL
#define BADPARAMETER_ERR    PASTIX_ERR_BADPARAMETER
#define FILE_ERR            PASTIX_ERR_FILE
#define BAD_DEFINE_ERR      PASTIX_ERR_BAD_DEFINE
#define INTEGER_TYPE_ERR    PASTIX_ERR_INTEGER_TYPE
#define IO_ERR              PASTIX_ERR_IO
#define MATRIX_ERR          PASTIX_ERR_MATRIX
#define FLOAT_TYPE_ERR      PASTIX_ERR_FLOAT_TYPE
#define STEP_ORDER_ERR      PASTIX_ERR_STEP_ORDER
#define MPI_ERR             PASTIX_ERR_MPI


#endif
