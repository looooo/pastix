/**
 *
 * @file cost.h
 *
 * PaStiX analyse headers for the cost matrix arrays functions.
 *
 * @copyright 1998-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @date 2024-07-05
 *
 * @addtogroup blend_dev_cost
 * @{
 *    This module contains all subroutines to initialize the cost arrays for a
 *    single matrix that will be used in the proportionnal mapping algorithm, as
 *    well as the simulation of the numerical factorization that defines the
 *    final mapping.
 *
 **/
#ifndef _cost_h_
#define _cost_h_

/**
 * @brief Arrays of double to store the cost of each element in the matrix
 */
typedef struct cost_matrix_s {
    double *blokcost; /**< Cost of the update generated by this block
                           for off-diagonal block, fact+solve otherwise */
    double *cblkcost; /**< Cost of all the operations linked to a panel */
} CostMatrix;

void        costMatrixInit ( CostMatrix *costmtx );
void        costMatrixExit ( CostMatrix *costmtx );
CostMatrix *costMatrixBuild( const symbol_matrix_t *symbmtx,
                             pastix_coeftype_t      flttype,
                             pastix_factotype_t     factotype );

#endif /* _cost_h_ */

/**
 * @}
 */
