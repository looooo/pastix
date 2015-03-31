/**
 *
 * @file symbol_cost.h
 *
 *  PaStiX symbol structure routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author David Goudin
 * @author Francois Pelegrin
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @date 2013-06-24
 *
 **/
#ifndef _SYMBOL_COST_H_
#define _SYMBOL_COST_H_

typedef struct symbol_function_s {
    double (*diag     )(pastix_int_t);
    double (*trsm     )(pastix_int_t, pastix_int_t);
    double (*update   )(pastix_int_t, pastix_int_t);
    double (*blkupdate)(pastix_int_t, pastix_int_t, pastix_int_t);
} symbol_function_t;

extern symbol_function_t flopstable[2][4];
extern symbol_function_t perfstable[2][4];

#endif /* _SYMBOL_COST_H_ */
