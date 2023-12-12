/**
 *
 * @file pastix_rhs.c
 *
 * @copyright 2004-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Vincent Bridonneau
 * @author Alycia Lisito
 * @date 2023-07-21
 *
 **/
#include "common.h"
#include "bcsc/bvec.h"
#include <lapacke.h>

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Initialize an RHS data structure.
 *
 *******************************************************************************
 *
 * @param[inout] B_ptr
 *          On entry, an allocated pastix_rhs_t data structure.
 *          On exit, the data is initialized to be used by the pastix_subtask_*
 *          functions.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixRhsInit( pastix_rhs_t  *B_ptr )
{
    pastix_rhs_t B;

    if ( B_ptr == NULL ) {
        pastix_print_error( "pastixRhsInit: wrong B parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }

    *B_ptr = malloc( sizeof(struct pastix_rhs_s) );
    B = *B_ptr;

    B->allocated = -1;
    B->flttype   = PastixPattern;
    B->m         = -1;
    B->n         = -1;
    B->ld        = -1;
    B->b         = NULL;
    B->cblkb     = NULL;
    B->rhs_comm  = NULL;
    B->Ploc2Pglob  = NULL;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Cleanup an RHS data structure.
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
pastixRhsFinalize( pastix_rhs_t B )
{
    if ( B == NULL ) {
        pastix_print_error( "pastixRhsFinalize: wrong B parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( B->b != NULL ) {
        if ( B->allocated > 0 ) {
            free( B->b );
        }
        else {
            pastix_print_warning( "Calling pastixRhsFinalize before restoring the ordering of vector b.\n"
                                  "Please call:\n"
                                  "  pastix_subtask_applyorder( pastix_data, flttype, PastixDirBackward, m, n,\n"
                                  "                             b, ldb, Bp );\n"
                                  "prior to this call to restore it.\n" );
        }
    }

    if ( B->cblkb != NULL ) {
        memFree_null( B->cblkb );
    }

    if ( B->Ploc2Pglob != NULL ) {
        memFree_null( B->Ploc2Pglob );
    }
    if ( B->rhs_comm != NULL ) {
        memFree_null( B->rhs_comm );
    }
    free( B );
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Reduces the precision of an RHS.
 *
 *******************************************************************************
 *
 * @param[in] dB
 *          The allocated pastix_rhs_t data structure to convert to lower
 *          precision.
 *
 * @param[out] sB
 *          On entry, an allocated pastix_rhs_t data structure.
 *          On exit, the reduced precision pastix_rhs_t of dB.
 *          If sB->allocated == -1 on entry, the internal b vector is
 *          automatically allocated by the function.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixRhsDoubletoSingle( const pastix_rhs_t dB,
                         pastix_rhs_t       sB )
{
    int rc;
    int tofree = 0;

    /* Generates halved-precision vector */
    if ( ( dB->flttype != PastixComplex64 ) &&
         ( dB->flttype != PastixDouble    ) )
    {
        pastix_print_error( "bvecDoubletoSingle: Invalid float type for mixed-precision" );
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( sB->allocated == -1 ) {
        size_t size = dB->ld * dB->n;

        memcpy( sB, dB, sizeof( struct pastix_rhs_s ) );

        sB->allocated = 1;
        sB->flttype   = dB->flttype - 1;
        sB->b         = malloc( size * pastix_size_of( sB->flttype ) );
        sB->rhs_comm  = NULL;
        tofree        = 1;
    }
    assert( sB->allocated >= 0 );
    assert( sB->flttype == (dB->flttype - 1) );
    assert( sB->b       != NULL );
    assert( sB->m       == dB->m );
    assert( sB->n       == dB->n );

    switch( dB->flttype ) {
    case PastixComplex64:
        rc = LAPACKE_zlag2c_work( LAPACK_COL_MAJOR, dB->m, dB->n, dB->b, dB->ld, sB->b, sB->ld );
        break;
    case PastixDouble:
        rc = LAPACKE_dlag2s_work( LAPACK_COL_MAJOR, dB->m, dB->n, dB->b, dB->ld, sB->b, sB->ld );
        break;
    default:
        rc = 1;
        pastix_print_error( "bvecDoubletoSingle: Invalid input float type for mixed-precision" );
    }

    if ( rc ) {
        if ( tofree ) {
            free( dB->b );
            dB->b = NULL;
        }
        return PASTIX_ERR_INTERNAL;
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Increases the precision of an RHS.
 *
 *******************************************************************************
 *
 * @param[in] sB
 *          The allocated pastix_rhs_t data structure to convert to higher
 *          precision.
 *
 * @param[out] dB
 *          On entry, an allocated pastix_rhs_t data structure.
 *          On exit, the increased precision pastix_rhs_t of sB.
 *          If dB->allocated == -1 on entry, the internal b vector is
 *          automatically allocated by the function.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixRhsSingleToDouble( const pastix_rhs_t sB,
                         pastix_rhs_t       dB )
{
    int rc;
    int tofree = 0;

    /* Frees halved-precision vector */
    if ( ( sB->flttype != PastixComplex32 ) &&
         ( sB->flttype != PastixFloat     ) )
    {
        pastix_print_error( "bvecSingleToDouble: Invalid input float type for mixed-precision" );
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( dB->allocated == -1 ) {
        size_t size = sB->ld * sB->n;

        memcpy( dB, sB, sizeof( struct pastix_rhs_s ) );

        dB->allocated = 1;
        dB->flttype   = sB->flttype + 1;
        dB->b         = malloc( size * pastix_size_of( dB->flttype ) );
        dB->rhs_comm  = NULL;
        tofree        = 1;
    }
    assert( dB->allocated >= 0 );
    assert( dB->flttype == (sB->flttype + 1) );
    assert( dB->b       != NULL );
    assert( dB->m       == sB->m );
    assert( dB->n       == sB->n );

    switch( sB->flttype ) {
    case PastixComplex32:
        rc = LAPACKE_clag2z_work( LAPACK_COL_MAJOR, sB->m, sB->n, sB->b, sB->ld, dB->b, dB->ld );
        break;
    case PastixFloat:
        rc = LAPACKE_slag2d_work( LAPACK_COL_MAJOR, sB->m, sB->n, sB->b, sB->ld, dB->b, dB->ld );
        break;
    default:
        rc = 1;
        pastix_print_error( "bvecSingleToDouble: Invalid float type for mixed-precision" );
    }

    if ( rc ) {
        if ( tofree ) {
            free( sB->b );
            sB->b = NULL;
        }
        return PASTIX_ERR_INTERNAL;
    }

    return PASTIX_SUCCESS;
}
