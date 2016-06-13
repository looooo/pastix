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

int  (*coeftabDiff[4])(const SolverMatrix*, SolverMatrix*);
void (*coeftabUncompress[4])(SolverMatrix*);
void (*coeftabCompress[4])(SolverMatrix*);

#endif /* _COEFTAB_H_ */
