/**
 * @file coeftab.h
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author David Goudin
 * @author Pascal Henon
 * @author Francois Pellegrini
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#ifndef _COEFTAB_H_
#define _COEFTAB_H_

#include "sopalin/coeftab_z.h"
#include "sopalin/coeftab_c.h"
#include "sopalin/coeftab_d.h"
#include "sopalin/coeftab_s.h"

void
coeftabInit( const pastix_data_t *pastix_data,
             int fakefillin, int factoLU );

void
coeftabExit( SolverMatrix *solvmtx );

typedef pastix_int_t (*coeftab_fct_compress_t)  ( SolverMatrix * );
typedef void         (*coeftab_fct_uncompress_t)( SolverMatrix * );
typedef pastix_int_t (*coeftab_fct_memory_t)    ( const SolverMatrix * );
typedef int          (*coeftab_fct_diff_t)      ( const SolverMatrix *, SolverMatrix * );

coeftab_fct_diff_t       coeftabDiff[4];
coeftab_fct_memory_t     coeftabMemory[4];
coeftab_fct_compress_t   coeftabCompress[4];
coeftab_fct_uncompress_t coeftabUncompress[4];

#endif /* _COEFTAB_H_ */
