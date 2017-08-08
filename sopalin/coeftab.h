/**
 *
 * @file coeftab.h
 *
 * PaStiX coefficient array routines header.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _COEFTAB_H_
#define _COEFTAB_H_

#include "sopalin/coeftab_z.h"
#include "sopalin/coeftab_c.h"
#include "sopalin/coeftab_d.h"
#include "sopalin/coeftab_s.h"

void coeftabInit( pastix_data_t *pastix_data,
                  int factoLU );
void coeftabExit( SolverMatrix *solvmtx );

typedef pastix_int_t (*coeftab_fct_compress_t)  ( SolverMatrix * );
typedef void         (*coeftab_fct_uncompress_t)( SolverMatrix * );
typedef pastix_int_t (*coeftab_fct_memory_t)    ( const SolverMatrix * );
typedef int          (*coeftab_fct_diff_t)      ( const SolverMatrix *, SolverMatrix * );

coeftab_fct_diff_t       coeftabDiff[4];
coeftab_fct_memory_t     coeftabMemory[4];
coeftab_fct_compress_t   coeftabCompress[4];
coeftab_fct_uncompress_t coeftabUncompress[4];

#endif /* _COEFTAB_H_ */
