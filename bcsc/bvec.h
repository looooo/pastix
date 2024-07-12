/**
 *
 * @file bvec.h
 *
 * @copyright 2004-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Vincent Bridonneau
 * @author Alycia Lisito
 * @date 2024-07-05
 *
 * @addtogroup bcsc
 * @{
 *
 **/
#ifndef _bvec_h_
#define _bvec_h_

#include <stdlib.h>

/**
 * @brief Tags used in MPI communications.
 */
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
    pastix_int_t idxcnt; /**< Amount of indexes of b or x which will be exchanged. */
    pastix_int_t valcnt; /**< Amount of values of b or x which will be exchanged.  */
} bvec_data_amount_t;

/**
 * @brief Informations of the data exchanged with other processes.
 */
typedef struct bvec_proc_comm_s
{
    bvec_data_amount_t  nsends;      /**< Number of indexes and values to send to clustnum.          */
    bvec_data_amount_t  nrecvs;      /**< Number of indexes and values to receive from clustnum.     */
    pastix_int_t       *send_idxbuf; /**< Array of indexes of b or x to exchange.                    */
    void               *send_valbuf; /**< Array of values of b or x to exchange.                     */
                                     /*   The indexes and values array are sorted the following way: */
                                     /*   k is incremented -> k += dofk * j.                         */
                                     /*   indexes[k] = i.                                            */
                                     /*   array[k]                -> corresponds to (i,     0).      */
                                     /*   array[k + dofk * j]     -> corresponds to (i,     j).      */
                                     /*   array[k + h + dofk * j] -> corresponds to (i + h, j).      */
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

#if defined( PASTIX_WITH_MPI )
int bvec_handle_comm_init( const pastix_data_t *pastix_data,
                           pastix_rhs_t         Pb );
int bvec_handle_comm_exit( bvec_handle_comm_t *rhs_comm );

pastix_int_t bvec_glob2Ploc( const pastix_data_t *pastix_data,
                             pastix_int_t         ig );
pastix_int_t bvec_Pglob2loc( const pastix_data_t *pastix_data,
                             const pastix_int_t  *glob2loc,
                             pastix_int_t         igp );
int bvec_compute_Ploc2Pglob( pastix_data_t *pastix_data,
                             pastix_rhs_t   Pb );

int bvec_exchange_amount_rep( bvec_handle_comm_t *rhs_comm );
int bvec_exchange_amount_dst( pastix_data_t *pastix_data,
                              pastix_dir_t   dir,
                              pastix_int_t   m,
                              pastix_int_t   nrhs,
                              pastix_rhs_t   Pb );
#endif

#endif /* _bvec_h_ */
/**
 * @}
 */
