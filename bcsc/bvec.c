/**
 *
 * @file bvec.c
 *
 * @copyright 2004-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Vincent Bridonneau
 * @author Alycia Lisito
 * @date 2022-12-05
 *
 **/
#include "common.h"
#include "bvec.h"

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Allocate a vector
 *
 *******************************************************************************
 *
 * @param[in] size
 *          The size of the vector
 *
 *******************************************************************************
 *
 * @return The allocated vector
 *
 *******************************************************************************/
void *bvec_malloc( size_t size )
{
    void *x = NULL;
    MALLOC_INTERN(x, size, char);
    return x;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Free a vector
 *
 *******************************************************************************
 *
 * @param[inout] x
 *          The vector to be free
 *
 *******************************************************************************/
void bvec_free( void *x )
{
    memFree_null(x);
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Cleanup the communicator part of the bvec data structure.
 *
 *******************************************************************************
 *
 * @param[inout] B
 *          On entry, the initialized pastix_rhs_t data structure.
 *          On exit, the structure is destroyed and should no longer be used.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
bvecHandleCommExit( bvec_handle_comm_t *rhs_comm )
{
    int c;
    int clustnbr = rhs_comm->clustnbr;

    for ( c = 0; c < clustnbr; c++ ) {
        if ( rhs_comm->data_comm[c].send.idxbuf != NULL ) {
            free( rhs_comm->data_comm[c].send.idxbuf );
        }
        if ( rhs_comm->data_comm[c].send.valbuf != NULL ) {
            free( rhs_comm->data_comm[c].send.valbuf );
        }
        if ( rhs_comm->data_comm[c].recv.idxbuf != NULL ) {
            free( rhs_comm->data_comm[c].recv.idxbuf );
        }
        if ( rhs_comm->data_comm[c].recv.valbuf != NULL ) {
            free( rhs_comm->data_comm[c].recv.valbuf );
        }
    }
    return PASTIX_SUCCESS;
}
