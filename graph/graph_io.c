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
    pastix_spm_t spm;
    FILE *stream;

    assert( pastix_data->procnbr == 1 );

    PASTIX_FOPEN(stream, "graphname","r");
    spmLoad( &spm, stream );
    fclose(stream);

    spmConvert( PastixCSC, &spm );
    graph->gN       = spm.gN;
    graph->n        = spm.n;
    graph->dof      = spm.dof;
    assert( spm.dof == 1 );
    graph->colptr   = spm.colptr;
    graph->rows     = spm.rowptr;
    graph->loc2glob = spm.loc2glob;

    (void)pastix_data;
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
    pastix_spm_t spm;
    FILE *stream;

    assert( pastix_data->procnbr == 1 );
    spm.n   = graph->n;
    spm.nnz = graph->colptr[ graph->n ] - graph->colptr[ 0 ];
    spm.dof = graph->dof;
    assert( spm.dof == 1 );
    spm.colptr   = graph->colptr;
    spm.rowptr   = graph->rows;
    spm.loc2glob = graph->loc2glob;
    spm.dofs     = NULL;

    spmUpdateFields( &spm );

    PASTIX_FOPEN(stream, "graphgen","r");
    spmSave( &spm, stream );
    fclose(stream);

    (void)pastix_data;
}
