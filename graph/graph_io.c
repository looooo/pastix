/**
 *
 * @file graph_io.c
 *
 *  PaStiX graph routines
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
#include "common.h"
#include "graph.h"
#include "spm.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * graphLoad - Load a graph from a file
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pointer to the solver instance to get options as rank,
 *          communicators, ...
 *
 * @param[in,out] graph
 *          The graph structure to store the loaded graph.
 *
 *******************************************************************************/
void graphLoad( const pastix_data_t  *pastix_data,
                      pastix_graph_t *graph )
{
    pastix_int_t  procnum  = pastix_data->procnum;
    pastix_int_t  ncol = 0;
    pastix_int_t *colptr, *rows;
    int dof = 1;

    if (procnum == 0) {
        FILE *stream;
        //TODO
        PASTIX_FOPEN(stream, "graphname","r");
        csc_load( &ncol, &colptr, &rows, NULL, NULL, &dof, stream);
        fclose(stream);
    }

    /* Number of columns */
    MPI_Bcast(&ncol, 1, PASTIX_MPI_INT, 0, pastix_data->pastix_comm);

        /* Colptr */
    if (procnum != 0) {
        MALLOC_INTERN(colptr, ncol+1, pastix_int_t);
    }
    MPI_Bcast(colptr, ncol+1, PASTIX_MPI_INT,
              0, pastix_data->pastix_comm);

    /* Rows */
    if  (procnum != 0) {
        MALLOC_INTERN(rows, colptr[ncol]-colptr[0], pastix_int_t);
    }
    MPI_Bcast(rows, colptr[ncol]-colptr[0], PASTIX_MPI_INT,
              0, pastix_data->pastix_comm);

    graph->gN       = ncol;
    graph->n        = ncol;
    graph->dof      = dof;
    graph->colptr   = colptr;
    graph->rows     = rows;
    graph->loc2glob = NULL;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * graphSave - Save a graph to file
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pointer to the solver instance to get options as rank,
 *          communicators, ...
 *
 * @param[in] graph
 *          The graph structure to store the loaded graph.
 *
 *******************************************************************************/
void graphSave( const pastix_data_t  *pastix_data,
                const pastix_graph_t *graph )
{
    pastix_int_t procnum = pastix_data->procnum;
    FILE *stream;

    if (procnum == 0)
    {
        PASTIX_FOPEN(stream, "cscgen", "w");
        // TODO
        csc_save( graph->n,
                  graph->colptr,
                  graph->rows,
                  0, NULL,
                  graph->dof, stream );
        fclose(stream);
    }
}
