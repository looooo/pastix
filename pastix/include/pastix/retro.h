#ifndef PASTIX_RETRO_H
#define PASTIX_RETRO_H

#define PASTIX_INT     pastix_int_t
#define pastix_float_t double

#define ASSERTDBG( test, mod ) assert(test)
#define ASSERT( test, mod ) assert(test)


/* #define PASTIX_PREFIX_F(x)   _PASTIX_ ## x */
/* #define PASTIX_EXTERN_F(x)   x */
/* #define PASTIX_PREFIX(x)     _PASTIX_ ## x */
/* #define PASTIX_EXTERN(x)     x */

#ifdef PASTIX_HAVE_MPI
#  define CHECK_MPI(call) do {                                          \
        int error_code;                                                 \
        error_code = call;                                              \
        if (error_code != MPI_SUCCESS) {                                \
                                                                        \
            char error_string[MPI_MAX_ERROR_STRING];                    \
            int length_of_error_string, error_class;                    \
            int my_rank = -1;                                           \
                                                                        \
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);                    \
            MPI_Error_class(error_code, &error_class);                  \
            MPI_Error_string(error_class, error_string,                 \
                             &length_of_error_string);                  \
            fprintf(stderr, "%3d: %s\n", my_rank, error_string);        \
            MPI_Error_string(error_code, error_string,                  \
                             &length_of_error_string);                  \
            fprintf(stderr, "%3d: %s\n", my_rank, error_string);        \
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                         \
        }                                                               \
    } while(0)
#else
#  define CHECK_MPI(call) call
#endif
#if (defined PASTIX_HAVE_MPI && defined PASTIX_HAVE_PTHREAD)
#    define CHECK_THREAD_LEVEL(THREAD_MODE)                             \
  do {                                                                  \
    int      provided;                                                  \
                                                                        \
    CHECK_MPI(MPI_Query_thread(&provided));                             \
    if (THREAD_MODE == API_THREAD_FUNNELED)                             \
      {                                                                 \
        switch(provided) {                                              \
        case MPI_THREAD_SINGLE:                                         \
          errorPrint("This run only supports MPI_THREAD_SINGLE\n"       \
                     "  either use -DFORCE_NOSMP,\n"                    \
                     "  change your MPI Library\n"                      \
                     "  or check that MPI_Init_thread"                  \
                     " is correctly called\n");                         \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        case MPI_THREAD_FUNNELED:                                       \
        case MPI_THREAD_SERIALIZED:                                     \
        case MPI_THREAD_MULTIPLE:                                       \
          break;                                                        \
        default:                                                        \
          errorPrint("provided thread level support is unknown");       \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        }                                                               \
      }                                                                 \
    else                                                                \
      {                                                                 \
        switch(provided) {                                              \
        case MPI_THREAD_SINGLE:                                         \
          errorPrint("This run only supports MPI_THREAD_SINGLE\n"       \
                     "  either use -DFORCE_NOSMP,\n"                    \
                     "  change your MPI Library\n"                      \
                     "  or check that MPI_Init_thread"                  \
                     " is correctly called\n");                         \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        case MPI_THREAD_FUNNELED:                                       \
          errorPrint("This run only supports MPI_THREAD_FUNNELED\n"     \
                     "  either use API_THREAD_FUNNELED,\n"              \
                     "  change your MPI Library\n"                      \
                     "  or check that MPI_Init_thread"                  \
                     " is correctly called\n");                         \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        case MPI_THREAD_SERIALIZED:                                     \
          errorPrint("This run only supports MPI_THREAD_SERIALIZED\n"   \
                     "  either use API_THREAD_FUNNELED,\n"              \
                     "  change your MPI Library\n"                      \
                     "  or check that MPI_Init_thread"                  \
                     " is correctly called\n");                         \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        case MPI_THREAD_MULTIPLE:                                       \
          break;                                                        \
        default:                                                        \
          errorPrint("provided thread level support is unknown");       \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        }                                                               \
      }                                                                 \
  } while(0)
#else
#  define CHECK_THREAD_LEVEL(THREAD_MODE) do {  \
    THREAD_MODE = THREAD_MODE;                  \
  } while (0)
#endif
#define CONJ_FLOAT(val) (val)
#define ABS_FLOAT(val)  fabs(val)

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
