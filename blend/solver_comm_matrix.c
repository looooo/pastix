/**
 *
 * @file solver_com_matrix.c
 *
 * PaStiX communication matrix handler.
 *
 * @copyright 1998-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Nolan Bredel
 * @date 2021-07-02
 *
 **/
#include "common/common.h"
#include "blend/solver_comm_matrix.h"

#if !defined(PASTIX_COMMUNICATION_MATRIX)
#error "PASTIX_COMMUNICATION_MATRIX option is disabled"
#endif

/**
 *******************************************************************************
 *
 * @brief Initialize the communication matrix.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver structure to which to attach the communication vector.
 *
 *******************************************************************************/
void
solverComMatrixInit( SolverMatrix *solvmtx )
{
    solvmtx->com_vector = calloc( solvmtx->clustnbr, sizeof( size_t ) );
}

/**
 *******************************************************************************
 *
 * @brief Free the communication matrix.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver structure to which to free the communication vector.
 *
 *******************************************************************************/
void
solverComMatrixExit( SolverMatrix *solvmtx )
{
    free( solvmtx->com_vector );
    solvmtx->com_vector = NULL;
}

/**
 *******************************************************************************
 *
 * @brief Gather the volume of communication and save it as a csv.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver structure with the communication vector to gather.
 *
 *******************************************************************************/
void
solverComMatrixGather( SolverMatrix *solvmtx )
{
    size_t *com_matrix = NULL;
    int     src;
    FILE *  file;

    MPI_Comm_rank( solvmtx->solv_comm, &src );
    if ( src == 0 ) {
        com_matrix = malloc( solvmtx->clustnbr * solvmtx->clustnbr * sizeof( size_t ) );
    }

    MPI_Gather( solvmtx->com_vector, solvmtx->clustnbr, MPI_UNSIGNED_LONG,
                com_matrix,          solvmtx->clustnbr, MPI_UNSIGNED_LONG,
                0, solvmtx->solv_comm );

    if ( src == 0 )
    {
        file = pastix_fopenw( ".", "com_matrix.csv", "w" );

        /* header */
        int p;
        for ( p = 0; p < solvmtx->clustnbr; p++ ) {
            fprintf( file, "R%d, ", p );
        }
        fprintf( file, "\n" );

        int i;
        for ( i = 0; i < ( solvmtx->clustnbr * solvmtx->clustnbr ); i++ ) {
            fprintf( file, "%lu, ", com_matrix[i] );
            if ( ( ( i + 1 ) % solvmtx->clustnbr ) == 0 ) {
                fprintf( file, "\n" );
            }
        }
        fclose(file);
    }
    free( com_matrix );
}
