/**
 *
 * @file queue.h
 *
 * PaStiX queue structure header.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @addtogroup blend_dev_queue
 * @{
 *    This module describes the queue structure used in the analyze part of the solver.
 *    The sorting is based on a balanced tree that is partially updated at
 *    insertion and suppression.
 *
 **/
#ifndef QUEUE_H
#define QUEUE_H

/**
 * @brief Queue item structure.
 */
typedef struct pastix_queue_item_s {
    double       key1;   /**< Key 1 of the element   */
    double       key2;   /**< Key 2 of the element   */
    pastix_int_t eltptr; /**< Pointer to the element */
} pastix_queue_item_t;

/**
 * @brief Queue structure.
 */
typedef struct pastix_queue_s {
    pastix_int_t         size;   /**< Allocated memory size          */
    pastix_int_t         used;   /**< Number of element in the queue */
    pastix_queue_item_t *elttab; /**< Array of the element           */
} pastix_queue_t;

int          pqueueInit(        pastix_queue_t *, pastix_int_t );
void         pqueueExit(        pastix_queue_t * );
pastix_int_t pqueueSize(  const pastix_queue_t * );
void         pqueueClear(       pastix_queue_t * );
void         pqueuePush2(       pastix_queue_t *, pastix_int_t, double, double );
pastix_int_t pqueueRead ( const pastix_queue_t * );
pastix_int_t pqueuePop2 (       pastix_queue_t *, double *, double * );
void         pqueuePrint( const pastix_queue_t * );

static inline void
pqueuePush1(pastix_queue_t *q, pastix_int_t elt, double key1) {
    pqueuePush2( q, elt, key1, 0. );
}

static inline pastix_int_t
pqueuePop(pastix_queue_t *q){
    return pqueuePop2(q, NULL, NULL);
}

static inline pastix_int_t
pqueuePop1(pastix_queue_t *q, double *key1){
    return pqueuePop2(q, key1, NULL);
}

#endif

/**
 * @}
 */
