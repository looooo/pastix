/**
 *
 * @file symbol_reorder.h
 *
 * PaStiX symbol structure routines
 *
 * @copyright 2004-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Vincent Bridonneau
 * @author Mathieu Faverge
 * @date 2024-07-05
 *
 * @addtogroup symbol_dev_reordering
 * @{
 *   @brief Functions to reorder the cblks of a matrix
 *
 **/
#ifndef _symbol_reorder_h_
#define _symbol_reorder_h_

void
symbol_reorder_cblk( const symbol_matrix_t *symbptr,
                     const symbol_cblk_t   *cblk,
                     pastix_order_t        *order,
                     const pastix_int_t    *levels,
                     pastix_int_t          *depthweight,
                     pastix_int_t           depthmax,
                     pastix_int_t           split_level,
                     pastix_int_t           stop_criterion );

void
symbol_reorder( pastix_data_t *pastix_data,
                pastix_int_t   maxdepth,
                pastix_int_t  *levels );

#endif /* _symbol_reorder_h_ */

/**
 * @}
 */
