/**
 *
 * @file solver_comm_matrix.h
 *
 * PaStiX communication matrix handler.
 *
 * @copyright 2021-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Nolan Bredel
 * @date 2021-07-02
 *
 **/
#ifndef _solver_comm_matrix_h_
#define _solver_comm_matrix_h_

#include "common/common.h"
#include "blend/solver.h"

#if defined(PASTIX_COMMUNICATION_MATRIX)
void solverComMatrixInit( SolverMatrix *solvmtx );
void solverComMatrixExit( SolverMatrix *solvmtx );
void solverComMatrixGather( SolverMatrix *solvmtx );

/**
 *******************************************************************************
 *
 * @brief Add the size of a communication
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] dest
 *          The destination node where the data goes.
 *
 * @param[in] size
 *          Size of the communication.
 *
 *******************************************************************************/
static inline void
solverCommMatrixAdd( SolverMatrix *solvmtx, int dest, size_t size )
{
    solvmtx->com_vector[dest] += size;
}

#else
static inline void
solverComMatrixInit( __attribute__((unused)) SolverMatrix *solvmtx ) {}

static inline void
solverComMatrixExit( __attribute__((unused)) SolverMatrix *solvmtx ) {}

static inline void
solverComMatrixGather( __attribute__((unused)) SolverMatrix *solvmtx ) {}

static inline void
solverCommMatrixAdd( __attribute__((unused)) SolverMatrix *solvmtx,
                     __attribute__((unused)) int           dest,
                     __attribute__((unused)) size_t        size )
{}
#endif

#endif /* _solver_comm_matrix_h_ */
