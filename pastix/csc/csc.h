/**
 *
 * @file csc.h
 *
 *  PaStiX csc routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _CSC_H_
#define _CSC_H_

/**
 * @ingroup pastix_csc
 * @struct pastix_csc_s - Csc structure.
 */
struct pastix_csc_s {
    pastix_int_t  gN;        /*< Global number of vertices                    */
    pastix_int_t  n;         /*< Number of local vertices                     */
    pastix_int_t *colptr;    /*< List of indirections to rows for each vertex */
    pastix_int_t *rows;      /*< List of edges for each vertex                */
    pastix_int_t *loc2glob;  /*< Corresponding numbering from local to global */
    int           ft;
    void         *avals;
};
typedef struct pastix_csc_s pastix_csc_t;


int
csc_load( pastix_int_t  *n,
          pastix_int_t **colptr,
          pastix_int_t **rows,
          int           *valtype,
          void         **values,
          int           *dof,
          FILE          *infile );

int
csc_save( pastix_int_t  n,
          pastix_int_t *colptr,
          pastix_int_t *rows,
          int           ft,
          void         *values,
          int           dof,
          FILE         *outfile );

int cscLoad( pastix_csc_t *csc, FILE *infile );
int cscSave( pastix_csc_t *csc, FILE *outfile );


#endif /* _CSC_H_ */
