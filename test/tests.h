/**
 *
 * @file tests.h
 *
 * Tests functions header.
 *
 * @copyright 2018-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#ifndef _tests_h_
#define _tests_h_

#include <stdio.h>

struct test_param_s {
    int n[3];       /**< Matrix size (min, max, step)            */
    int mode[3];    /**< Matrix generation mode (min, max, step) */
    int prank[3];   /**< Matrix rank percentage (min, max, step) */
    int method[3];  /**< Compression method (min, max, step)     */
    double tol_gen; /**< Tolerance for the matrix generation     */
    double tol_cmp; /**< Tolerance for the matrix compression    */
    double threshold; /**< Tolerance for the matrix compression    */
    FILE  *output;
};
typedef struct test_param_s test_param_t;

void testGetOptions( int argc, char **argv, test_param_t *params, double eps );

#endif /* _tests_h_ */
