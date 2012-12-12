/**
 *
 * @file pzpack.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.4.6
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @precisions normal z -> c d s
 *
 **/

#include <stdlib.h>
#include <sys/types.h>
#include "common.h"

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 * plasma_pzpack pack all extra elements at the end of the matrix
 *
 *      +---------------+
 *      |               |
 *      |               |
 *      |     A11       |
 *      |               |
 *      |               |
 *      +---------------+
 *      |     A21       |
 *      +---------------+
 *
 * This matrix is initially stored as (example of Column Major, it's
 * the same for row major. We just consider the transpose matrix) :
 *  A11(:,0), A21(:,0), A11(:,1), A21(:,1), ... 
 *
 * On exit, it's stored as follow.
 *  A11(:,:), A12(:,:)
 *******************************************************************************
 *
 * @param[in] plasma
 *         Plasma context
 *
 * Parameters given to the function:
 *     m [in]     Number of rows in matrix A
 *     n [in]     Number of columns in matrix A
 *     A [in,out] Matrix A to pack. (see above for entry and exit format)
 *     m0 [in]    Number of rows of A21
 *
 ******************************************************************************/
void plasma_pzpack(plasma_context_t *plasma)
{
    PLASMA_Complex64_t *A, *W, *Wl;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int m, n, m0;
    int i, m1, size, rank, start, end, bs, mod;

    plasma_unpack_args_6(m, n, A, m0, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    /* Quick return */
    if ( n <= 1 )
      return;

    m1 = m - m0;

    size = PLASMA_SIZE;
    rank = PLASMA_RANK;

    mod   = (n-1) % size;
    bs    = (n-1) / size;
    start = rank * bs;
    if ( rank < mod ) {
        bs++;
    }
    start += min( mod, rank );

    W  = (PLASMA_Complex64_t*)plasma_private_alloc(plasma, (m0*bs), PlasmaComplexDouble);
    Wl = (PLASMA_Complex64_t*)plasma_private_alloc(plasma, m1,      PlasmaComplexDouble);

    /* Save leftover pieces that are otherwise going to be overwritten */
    CORE_zlacpy( PlasmaUpperLower, m0, bs, &(A[(int64_t)start*m+m1]), m, W, m0 );

    /* Pack A */
    end = ((n-1) / size) * size + 1;
    for(i=rank+1; i<end; i+=size) {
        memcpy( Wl, &(A[i*m]), m1*sizeof(PLASMA_Complex64_t));
        plasma_barrier(plasma);
        memcpy( &(A[i*m1]), Wl, m1*sizeof(PLASMA_Complex64_t));
    }

    if ( rank < (n - end)) {
        i = end + rank;
        memcpy( Wl, &(A[i*m]), m1*sizeof(PLASMA_Complex64_t));
        plasma_barrier(plasma);
        memcpy( &(A[i*m1]), Wl, m1*sizeof(PLASMA_Complex64_t));
    }
    else
        plasma_barrier(plasma);

    /* Restore leftover pieces */
    CORE_zlacpy( PlasmaUpperLower, m0, bs, W, m0, &(A[(int64_t)m1*n+start*m0]), m0 );

    plasma_private_free(plasma, W);
    plasma_private_free(plasma, Wl);
}


/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 * plasma_pzunpack unpack all extra elements from the end of the matrix
 *
 *      +---------------+
 *      |               |
 *      |               |
 *      |     A11       |
 *      |               |
 *      |               |
 *      +---------------+
 *      |     A21       |
 *      +---------------+
 *
 * This matrix is initially stored as (example of Column Major, it's
 * the same for row major. We just consider the transpose matrix) :
 *  A11(:,:), A12(:,:)
 *
 * On exit, it's stored as follow.
 *  A11(:,0), A21(:,0), A11(:,1), A21(:,1), ... 
 *
 *******************************************************************************
 *
 * @param[in] plasma
 *         Plasma context
 *
 * Parameters given to the function:
 *     m [in]     Number of rows in matrix A
 *     n [in]     Number of columns in matrix A
 *     A [in,out] Matrix A to unpack. (see above for entry and exit format)
 *     m0 [in]    Number of rows of A21
 *
 ******************************************************************************/
void plasma_pzunpack(plasma_context_t *plasma)
{
    PLASMA_Complex64_t *A, *W, *Wl;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int m, n, m0;
    int i, m1, size, rank, start, end, bs, mod;

    plasma_unpack_args_6(m, n, A, m0, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    /* Quick return */
    if ( n <= 1 )
      return;

    m1 = m - m0;

    size = PLASMA_SIZE;
    rank = PLASMA_RANK;

    mod   = (n-1) % size;
    bs    = (n-1) / size;
    start = rank * bs;
    if ( rank < mod ) {
        bs++;
    }
    start += min( mod, rank );

    W  = (PLASMA_Complex64_t*)plasma_private_alloc(plasma, (m0*bs), PlasmaComplexDouble);
    Wl = (PLASMA_Complex64_t*)plasma_private_alloc(plasma, m1,      PlasmaComplexDouble);

    /* Save leftover pieces that are otherwise going to be overwritten */
    CORE_zlacpy( PlasmaUpperLower, m0, bs, &(A[(int64_t)m1*n+start*m0]), m0, W, m0 );

    /* Unpack A */
    end = ((n-1) % size) ;
    for(i=n-1-rank; i>end; i-=size) {
        memcpy( Wl, &(A[i*m1]), m1*sizeof(PLASMA_Complex64_t));
        plasma_barrier(plasma);
        memcpy( &(A[i*m]), Wl, m1*sizeof(PLASMA_Complex64_t));
    }

    if ( rank < end ) {
        i = rank+1;
        memcpy( Wl, &(A[i*m1]), m1*sizeof(PLASMA_Complex64_t));
        plasma_barrier(plasma);
        memcpy( &(A[i*m]), Wl, m1*sizeof(PLASMA_Complex64_t));
    }
    else
        plasma_barrier(plasma);

    /* Restore leftover pieces */
    CORE_zlacpy( PlasmaUpperLower, m0, bs, W, m0, &(A[(int64_t)start*m+m1]), m );

    plasma_private_free(plasma, W);
    plasma_private_free(plasma, Wl);
}
