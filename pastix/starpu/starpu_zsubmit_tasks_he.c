/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
/*
 * File: starpu_zsubmit_tasks.c
 *
 * Functions computing operations for reffinement methods
 *
 */

#ifdef CHOL_SOPALIN
#undef CHOL_SOPALIN
#endif
#define HERMITIAN
#include "starpu_zsubmit_tasks.c"
