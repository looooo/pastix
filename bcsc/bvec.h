/**
 *
 * @file bvec.h
 *
 * @copyright 2004-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Vincent Bridonneau
 * @author Alycia Lisito
 * @date 2022-12-05
 *
 * @addtogroup bcsc
 * @{
 *
 **/
#ifndef _bvec_h_
#define _bvec_h_

#include <stdlib.h>

typedef enum bvec_tag_ {
    PastixTagAmount,
    PastixTagIndexes,
    PastixTagValues,
} bvec_tag_e;

/**
 * @brief Information about the amount of data exchanged to permute the pivots.
 */
typedef struct bvec_data_amount_s
{
    /* WARNING: the two counters must remain at the head of the structure to exchange
       the volume of the communications between the processes. */
    pastix_int_t  idxcnt; /**< Amount of indexes of b or x which will be exchanged. */
    pastix_int_t  valcnt; /**< Amount of values of b or x which will be exchanged.  */
    pastix_int_t *idxbuf; /**< Array of indexes of b or x to exchange.              */
    void         *valbuf; /**< Array of values of b or x to exchange.
                               The indexes and values array are sorted the following way:
                               k is incremented -> k += dofk * j.
                               indexes[k] = i.
                               array[k]                -> corresponds to (i,     0).
                               array[k + dofk * j]     -> corresponds to (i,     j).
                               array[k + h + dofk * j] -> corresponds to (i + h, j).  */
} bvec_data_amount_t;

/**
 * @brief Informations of the data exchanged with other processes.
 */
typedef struct bvec_proc_comm_s
{
    bvec_data_amount_t send; /**< Number of indexes and values to send to clustnum.      */
    bvec_data_amount_t recv; /**< Number of indexes and values to receive from clustnum. */
} bvec_proc_comm_t;

/**
 * @brief Structure to manage communications with distributed rhs.
 */
typedef struct bvec_handle_comm_s
{
    pastix_int_t      clustnbr;     /**< Number of processes in the communicator.                                       */
    pastix_int_t      clustnum;     /**< ID of the current process in the communicator.                                 */
    PASTIX_Comm       comm;         /**< PaStiX MPI communicator used for the ordering step.                            */
    pastix_coeftype_t flttype;      /**< valtab datatype: PastixFloat, PastixDouble, PastixComplex32 or PastixComplex64 */
    pastix_int_t      max_idx;      /**< Maximum amount of indexes received, used to allocate the receiving buffer.     */
    pastix_int_t      max_val;      /**< Maximum amount of values received, used to allocate the receiving buffer.      */
    bvec_proc_comm_t  data_comm[1]; /**< Array of size clustnbr.                                                        */
                                    /* data_comm[c]: contains the data clustnum has to send to c in                     */
                                    /*               the distributed case and the amount of data                        */
                                    /*               clustnum will receive from c.                                      */
                                    /* data_comm[clustnum]: contains the data clustnum will send                        */
                                    /*                      to the other processors in the replicated                   */
                                    /*                      case and is NULL in the distributed case.                   */
} bvec_handle_comm_t;

void *bvec_malloc( size_t size );
void  bvec_free( void *x );
int bvec_handle_comm_init( const pastix_data_t *pastix_data,
                           pastix_rhs_t         Pb );
int bvec_handle_comm_exit( bvec_handle_comm_t *rhs_comm );

pastix_int_t bvec_glob2Ploc( const pastix_data_t *pastix_data,
                             pastix_int_t         ig );

int bvec_Ploc2Pglob( pastix_data_t *pastix_data,
                     pastix_rhs_t   Pb );

int bvec_exchange_amount_rep( bvec_handle_comm_t *rhs_comm );

#endif /* _bvec_h_ */
/**
 * @}
 */
