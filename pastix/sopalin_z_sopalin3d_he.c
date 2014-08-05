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
 File: sopalin3d.c

 sopalin 3d main program.

 Authors:
 Mathieu Faverge - faverge@labri.fr
 Xavier  Lacoste - lacoste@labri.fr
 Pierre  Ramet   - ramet@labri.fr

 Date:
 Version 0.0 - february 2003
 */

#define CHOL_SOPALIN
#ifdef SOPALIN_LU
#undef SOPALIN_LU
#endif
#include "sopalin3d.c"
