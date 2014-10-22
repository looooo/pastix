/**
 *
 * @file sopalin_data.h
 *
 *  PaStiX sopalin routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _SOPALIN_DATA_H_
#define _SOPALIN_DATA_H_

struct sopalin_data_s {
    SolverMatrix *solvmtx;
};
typedef struct sopalin_data_s sopalin_data_t;


void pastix_static_zpotrf( sopalin_data_t *sopalin_data );
void pastix_static_cpotrf( sopalin_data_t *sopalin_data );
void pastix_static_dpotrf( sopalin_data_t *sopalin_data );
void pastix_static_spotrf( sopalin_data_t *sopalin_data );

#endif /* _SOPALIN_DATA_H_ */


