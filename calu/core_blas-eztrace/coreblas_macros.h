/**
 *
 * @file coreblas_macros.h
 *
 *  PLASMA core_blas-eztrace module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#ifndef _COREBLAS_MACROS_H_
#define _COREBLAS_MACROS_H_

/* #define VERBOSE 1 */
#ifdef VERBOSE
#define FUNCTION_NAME fprintf(stderr, "Calling [%s]\n", __FUNCTION__)
#else
#define FUNCTION_NAME
#endif

/*
 * Macro for coreblas functions which returns void.
 * It generates an event with the BLAS/LAPACK name
 * when it starts and a STOP event when it's finished.
 */
#define FUNCTION_VOID(name, shname, rettype, args, params)	\
  rettype P##name args ;                                        \
  rettype name args {                                           \
    FUNCTION_NAME;                                              \
    EZTRACE_EVENT0(FUT_COREBLAS(shname));                       \
    P##name params;                                             \
    EZTRACE_EVENT0(FUT_COREBLAS(STOP));                         \
  }

/*
 * Macro for coreblas functions which returns something
 * different from void.
 * It generates an event with the BLAS/LAPACK name
 * when it starts and a STOP event when it's finished.
 */
#define FUNCTION_TYPE(name, shname, rettype, args, params)	\
  rettype P##name args ;                                        \
  rettype name args {                                           \
    rettype ret;                                                \
    FUNCTION_NAME;                                              \
    EZTRACE_EVENT0(FUT_COREBLAS(shname));                       \
    ret = P##name params;                                       \
    EZTRACE_EVENT0(FUT_COREBLAS(STOP));                         \
    return ret;                                                 \
  }

/*
 * Macro for quark wrapper used in coreblas library.
 * It generates an event with the BLAS/LAPACK name
 * of the function called, and an event to decrease
 * the number of tasks available.
 * At the end, the STOP event is generated.
 */
#define FUNCTION_QUARK(name, shname)                                    \
  void P##name (Quark *quark);                                          \
  void name(Quark *quark) {                                             \
    FUNCTION_NAME;                                                      \
    EZTRACE_EVENT2(FUT_COREBLAS(shname),                                \
                   QUARK_Get_Priority(quark),                           \
                   QUARK_Get_Sequence(quark));                          \
    P##name(quark);                                                     \
    EZTRACE_EVENT0(FUT_COREBLAS(STOP));                                 \
  }

#endif	/* _COREBLAS_MACROS_H_ */
