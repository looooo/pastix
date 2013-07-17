/**
 *
 * @file graph.h
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
#ifndef _GRAPH_H_
#define _GRAPH_H_

struct pastix_graph_s {
    pastix_int_t  gN;        /*< Global number of vertices                    */
    pastix_int_t  n;         /*< Number of local vertices                     */
    pastix_int_t *colptr;    /*< List of indirections to rows for each vertex */
    pastix_int_t *nnz;       /*< Number of edges for each vertex              */
    pastix_int_t *rows;      /*< List of edges for each vertex                */
    pastix_int_t *loc2glob;  /*< Corresponding numbering from local to global */
};

void graphLoad( const pastix_data_t  *pastix_data,
                pastix_graph_t *graph );
void graphSave( const pastix_data_t  *pastix_data,
                const pastix_graph_t *graph );

void graphBase(  pastix_graph_t *graph, int baseval );
void graphClean( pastix_graph_t *graph );

int  graphSymmetrize(       pastix_int_t    n,
                      const pastix_int_t   *ia,
                      const pastix_int_t   *ja,
                      const pastix_int_t   *loc2glob,
                            pastix_graph_t *newgraph );

int  graphPrepare(      pastix_data_t   *pastix_data,
                        pastix_int_t     n,
                  const pastix_int_t    *colptr,
                  const pastix_int_t    *rows,
                  const pastix_int_t    *loc2glob,
                        pastix_graph_t **graph );

int  graphIsolate(       pastix_int_t   n,
                   const pastix_int_t  *colptr,
                   const pastix_int_t  *rows,
                         pastix_int_t   isolate_n,
                         pastix_int_t  *isolate_list,
                         pastix_int_t **new_colptr,
                         pastix_int_t **new_rows,
                         pastix_int_t **new_perm,
                         pastix_int_t **new_invp );

int graphApplyPerm( const pastix_graph_t *graphA,
                    const pastix_int_t   *perm,
                          pastix_graph_t *graphPA );
#endif /* _GRAPH_H_ */
