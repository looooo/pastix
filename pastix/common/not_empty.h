#ifndef NOT_EMPTY_H
#define NOT_EMPTY_H
#include "redefine_functions.h"

#ifdef SOPALIN_LU
#ifndef CHOL_SOPALIN
#define CHOL_SOPALIN
#endif
#endif

#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU
#    define NOT_EMPTY(filename)                                \
  void PASTIX_PREFIX_F(ge ## _not_empty_ ## filename) (void){  \
    return;                                                    \
  }
#  else /* not SOPALIN_LU */
#    define NOT_EMPTY(filename)                                \
  void PASTIX_PREFIX_F(po ## _not_empty_ ## filename) (void){  \
    return;                                                    \
  }
#  endif /* not SOPALIN_LU */
#else /* not CHOL_SOPALIN */
#  ifdef HERMITIAN
#    define NOT_EMPTY(filename)                                \
  void PASTIX_PREFIX_F(he ## _not_empty_ ## filename) (void){  \
    return;                                                    \
  }
#  else /* not HERMITIAN */
#    define NOT_EMPTY(filename)                                \
  void PASTIX_PREFIX_F(sy ## _not_empty_ ## filename) (void){  \
    return;                                                    \
  }
#  endif /* not HERMITIAN */
#endif /* not CHOL_SOPALIN */
#endif /* not NOT_EMPTY_H */
