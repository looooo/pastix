/**
 *
 * @file zgetmi.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.5.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @precisions normal z -> c d s
 *
 **/

#include <stdlib.h>
#include <sys/types.h>
#include "common.h"
#include "zgecfi2.h"

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 *  PLASMA_zgecfi convert the matrice A in place from format f_in to
 *  format f_out
 *
 *******************************************************************************
 *
 * @param[in] m
 *         Number of rows of matrix A
 *
 * @param[in] n
 *         Number of columns of matrix A
 *
 * @param[in,out] A
 *         Matrix of size L*m*n
 *
 * @param[in] f_in
 *         Original format of the matrix A. Must be part of (PlasmaCM, PlasmaRM,
 *         PlasmaCCRB, PlasmaCRRB, PlasmaRCRB, PlasmaRRRB)
 *
 * @param[in] imb
 *         Number of rows of each block in original format
 *
 * @param[in] inb
 *         Number of columns of each block in original format
 *
 * @param[in] f_out
 *         Format requested for the matrix A. Must be part of (PlasmaCM, PlasmaRM,
 *         PlasmaCCRB, PlasmaCRRB, PlasmaRCRB, PlasmaRRRB)
 *
 * @param[in] omb
 *         Number of rows of each block in requested format
 *
 * @param[in] onb
 *         Number of columns of each block in requested format
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgecfi_Async
 *
 ******************************************************************************/
int PLASMA_zgecfi(int m, int n, PLASMA_Complex64_t *A,
                  PLASMA_enum f_in,  int imb, int inb,
                  PLASMA_enum f_out, int omb, int onb) 
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error(__func__, "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    plasma_sequence_create(plasma, &sequence);

    PLASMA_zgecfi_Async( m, n, A,
                         f_in,  imb, inb,
                         f_out, omb, onb,
                         sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);

    return status;
}


/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 *  PLASMA_zgecfi_Async convert the matrice A in place from format f_in to
 *  format f_out
 *
 *******************************************************************************
 *
 * @param[in] m
 *         Number of rows of matrix A
 *
 * @param[in] n
 *         Number of columns of matrix A
 *
 * @param[in,out] A
 *         Matrix of size L*m*n
 *
 * @param[in] f_in
 *         Original format of the matrix A. Must be part of (PlasmaCM, PlasmaRM,
 *         PlasmaCCRB, PlasmaCRRB, PlasmaRCRB, PlasmaRRRB)
 *
 * @param[in] imb
 *         Number of rows of each block in original format
 *
 * @param[in] inb
 *         Number of columns of each block in original format
 *
 * @param[in] f_out
 *         Format requested for the matrix A. Must be part of (PlasmaCM, PlasmaRM,
 *         PlasmaCCRB, PlasmaCRRB, PlasmaRCRB, PlasmaRRRB)
 *
 * @param[in] omb
 *         Number of rows of each block in requested format
 *
 * @param[in] onb
 *         Number of columns of each block in requested format
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgecfi
 *
 ******************************************************************************/
int PLASMA_zgecfi_Async(int m, int n, PLASMA_Complex64_t *A,
                        PLASMA_enum f_in,  int imb, int inb,
                        PLASMA_enum f_out, int omb, int onb,
                        PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_Complex64_t *W = NULL;
    int im1, in1, om1, on1;
    size_t A11, A21, A12, A22;

    /* Check Plasma context */
    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error(__func__, "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check arguments */
    if(  ( f_in  != PlasmaCM ) && ( f_in  != PlasmaRM ) 
         && ( f_in  != PlasmaCCRB ) && ( f_in  != PlasmaRRRB ) 
         && ( f_in  != PlasmaCRRB ) && ( f_in  != PlasmaRCRB ) )
    {
        plasma_error(__func__, "Input format unknown");
        return -4;
    }
    if(  ( f_out  != PlasmaCM ) && ( f_out  != PlasmaRM ) 
         && ( f_out  != PlasmaCCRB ) && ( f_out  != PlasmaRRRB ) 
         && ( f_out  != PlasmaCRRB ) && ( f_out  != PlasmaRCRB ) )
    {
        plasma_error(__func__, "Input format unknown");
        return -7;
    }

    /* quick return */
    if( (f_in == f_out) && ( (f_in == PlasmaCM) || (f_in == PlasmaRM))
        && (imb == omb) && ( inb == onb ) ) {
        return PLASMA_SUCCESS;
    }

    if ( (f_in == PlasmaCM) || (f_in == PlasmaRM) )
    {
        if ( (f_out == PlasmaCM) || (f_out == PlasmaRM) ){
            imb = omb = PLASMA_NB;
            inb = onb = PLASMA_NB;
        } else {
            imb = omb;
            inb = onb;
        }
    } 
    else if ( (f_out == PlasmaCM) || (f_out == PlasmaRM) ) 
    {
        omb = imb;
        onb = inb;
    }

    /* calculate number of full blocks */
    im1 = (m / imb) * imb;
    in1 = (n / inb) * inb;
    om1 = (m / omb) * omb;
    on1 = (n / onb) * onb;

    /* separate the four submatrices A11, A12, A21, A22 */
    if( f_in == PlasmaCM ) {
        if( om1 < m ) {
            plasma_static_call_6(plasma_pzpack,
                                 int,                 m,
                                 int,                 on1,
                                 PLASMA_Complex64_t*, A,
                                 int,                 (m-om1),
                                 PLASMA_sequence*,    sequence,
                                 PLASMA_request*,     request);
            if ( on1 < n) {
                plasma_static_call_6(plasma_pzpack,
                                     int,                 m,
                                     int,                 (n-on1),
                                     PLASMA_Complex64_t*, &(A[m*on1]),
                                     int,                 (m-om1),
                                     PLASMA_sequence*,    sequence,
                                     PLASMA_request*,     request);
            }
        }
    }
    else if ( f_in == PlasmaRM ) {
        if( on1 < n ) {
            plasma_static_call_6(plasma_pzpack,
                                 int,                 n,
                                 int,                 om1,
                                 PLASMA_Complex64_t*, A,
                                 int,                 (n-on1),
                                 PLASMA_sequence*,    sequence,
                                 PLASMA_request*,     request);
            if( om1 < m ) {
                plasma_static_call_6(plasma_pzpack,
                                     int,                 n,
                                     int,                 (m-om1),
                                     PLASMA_Complex64_t*, &(A[n*om1]),
                                     int,                 (n-on1),
                                     PLASMA_sequence*,    sequence,
                                     PLASMA_request*,     request);
            }
        }
    }

    /* blocked format to blocked format conversion with different block sizes */
    if( (f_in  != PlasmaCM) && (f_in  != PlasmaRM) &&
        (f_out != PlasmaCM) && (f_out != PlasmaRM) ) {
        if( (imb != omb) || (inb != onb) ) {
            if( (f_in == PlasmaRRRB)  || (f_out == PlasmaRRRB ) ) {
                PLASMA_zgecfi_Async(m, n, A, f_in, imb, inb, PlasmaRM,    1,   1,   sequence, request);
                PLASMA_zgecfi_Async(m, n, A, PlasmaRM,   1,   1,   f_out, omb, onb, sequence, request);
            }
            else {
                PLASMA_zgecfi_Async(m, n, A, f_in, imb, inb, PlasmaCM,    1,   1,   sequence, request);
                PLASMA_zgecfi_Async(m, n, A, PlasmaCM,   1,   1,   f_out, omb, onb, sequence, request);
            }
            return PLASMA_SUCCESS;
        }
    }

    if( (f_in == PlasmaCM) || (f_in == PlasmaCCRB) || (f_in == PlasmaCRRB) ) 
    {
        A11 = 0;
        A21 = im1*in1;
        A12 = m  *in1;
        A22 = m  *in1 + im1*(n-in1);
    }
    else
    {
        A11 = 0;
        A12 = im1*in1;
        A21 = im1*n;
        A22 = im1*n + in1*(m-im1);
    }

    switch ( f_in ) {
    case PlasmaCM :
        switch ( f_out ) {
        case PlasmaCM   : break;
        case PlasmaCCRB : ipt_call(cm2ccrb,   om1, on1, omb, onb); break;
        case PlasmaCRRB : ipt_call(cm2crrb,   om1, on1, omb, onb); break;
        case PlasmaRCRB : ipt_call(cm2rcrb,   om1, on1, omb, onb); break;
        case PlasmaRRRB : ipt_call(cm2rrrb,   om1, on1, omb, onb); break;
        case PlasmaRM   : ipt_call(cm2rm,     om1, on1, omb, onb); break;
        default: ;
        }
        break;
    case PlasmaCCRB:
        switch ( f_out ) {
        case PlasmaCM   : ipt_call(ccrb2cm,   im1, in1, imb, inb); break;
        case PlasmaCCRB : break;
        case PlasmaCRRB : ipt_cal2(ccrb2crrb, im1, in1, imb, inb); break;
        case PlasmaRCRB : ipt_call(ccrb2rcrb, im1, in1, imb, inb); break;
        case PlasmaRRRB : ipt_call(ccrb2rrrb, im1, in1, imb, inb); break;
        case PlasmaRM   : ipt_call(ccrb2rm,   im1, in1, imb, inb); break;
        default: ;
        }
        break;
    case PlasmaCRRB:
        switch ( f_out ) {
        case PlasmaCM   : ipt_call(crrb2cm,   im1, in1, imb, inb); break;
        case PlasmaCCRB : ipt_cal2(crrb2ccrb, im1, in1, imb, inb); break;
        case PlasmaCRRB : break;
        case PlasmaRCRB : ipt_call(crrb2rcrb, im1, in1, imb, inb); break;
        case PlasmaRRRB : ipt_call(crrb2rrrb, im1, in1, imb, inb); break;
        case PlasmaRM   : ipt_call(crrb2rm,   im1, in1, imb, inb); break;
        default: ;
        }
        break;
    case PlasmaRCRB:
        switch ( f_out ) {
        case PlasmaCM   : ipt_call(rcrb2cm,   im1, in1, imb, inb); break;
        case PlasmaCCRB : ipt_call(rcrb2ccrb, im1, in1, imb, inb); break;
        case PlasmaCRRB : ipt_call(rcrb2crrb, im1, in1, imb, inb); break;
        case PlasmaRCRB : break;
        case PlasmaRRRB : ipt_cal2(rcrb2rrrb, im1, in1, imb, inb); break;
        case PlasmaRM   : ipt_call(rcrb2rm,   im1, in1, imb, inb); break;
        default: ;
        }
        break;
    case PlasmaRRRB:
        switch ( f_out ) {
        case PlasmaCM   : ipt_call(rrrb2cm,   im1, in1, imb, inb); break;
        case PlasmaCCRB : ipt_call(rrrb2ccrb, im1, in1, imb, inb); break;
        case PlasmaCRRB : ipt_call(rrrb2crrb, im1, in1, imb, inb); break;
        case PlasmaRCRB : ipt_cal2(rrrb2rcrb, im1, in1, imb, inb); break;
        case PlasmaRRRB : break;
        case PlasmaRM   : ipt_call(rrrb2rm,   im1, in1, imb, inb); break;
        default: ;
        }
        break;
    case PlasmaRM:
        switch ( f_out ) {
        case PlasmaCM   : ipt_call(rm2cm,     om1, on1, omb, onb); break;
        case PlasmaCCRB : ipt_call(rm2ccrb,   om1, on1, omb, onb); break;
        case PlasmaCRRB : ipt_call(rm2crrb,   om1, on1, omb, onb); break;
        case PlasmaRCRB : ipt_call(rm2rcrb,   om1, on1, omb, onb); break;
        case PlasmaRRRB : ipt_call(rm2rrrb,   om1, on1, omb, onb); break;
        case PlasmaRM   : break;
        default: ;
        }
        break;
    default: ;
    }

    /* reorder block */
    if( (f_out == PlasmaCM) || (f_out == PlasmaCCRB) || (f_out == PlasmaCRRB) ) 
    {
        /* We need to swap A21 and A12 */
        if ( A21 > A12 ) {
            size_t sze1 = A21-A12;
            size_t sze2 = A22-A21;

            QUARK_Barrier(plasma->quark);
            //plasma_malloc(W, max( in1, on1), PLASMA_Complex64_t);
            W = (PLASMA_Complex64_t*)malloc( max( sze1, sze2 ) * sizeof(PLASMA_Complex64_t) );
            CORE_zswpab(0, sze1, sze2, &(A[A12]), W);
            free(W);
        }
    }
    else {
        /* We need to swap A21 and A12 */
        if ( A12 > A21 ) {
            size_t sze1 = A12-A21;
            size_t sze2 = A22-A12;

            QUARK_Barrier(plasma->quark);
            //plasma_malloc(W, max( in1, on1), PLASMA_Complex64_t);
            W = (PLASMA_Complex64_t*)malloc( max( sze1, sze2 ) * sizeof(PLASMA_Complex64_t) );
            CORE_zswpab(0, sze1, sze2, &(A[A21]), W);
            free(W);
        }
    }

    /* unseparate if output is not blocked */
    if( f_out == PlasmaCM ) {
        if( im1 < m ) {
            plasma_static_call_6(plasma_pzunpack,
                                 int,                 m,
                                 int,                 in1,
                                 PLASMA_Complex64_t*, A,
                                 int,                 (m-im1),
                                 PLASMA_sequence*,    sequence,
                                 PLASMA_request*,     request);
            if ( in1 < n) {
                plasma_static_call_6(plasma_pzunpack,
                                     int,                 m,
                                     int,                 (n-in1),
                                     PLASMA_Complex64_t*, &(A[m*in1]),
                                     int,                 (m-im1),
                                     PLASMA_sequence*,    sequence,
                                     PLASMA_request*,     request);
            }
        }
    }
    else if( f_out == PlasmaRM ) {
        if( in1 < n ) {
            plasma_static_call_6(plasma_pzunpack,
                                 int,                 n,
                                 int,                 im1,
                                 PLASMA_Complex64_t*, A,
                                 int,                 (n-in1),
                                 PLASMA_sequence*,    sequence,
                                 PLASMA_request*,     request);
            if( im1 < m ) {
                plasma_static_call_6(plasma_pzunpack,
                                     int,                 n,
                                     int,                 (m-im1),
                                     PLASMA_Complex64_t*, &(A[n*im1]),
                                     int,                 (n-in1),
                                     PLASMA_sequence*,    sequence,
                                     PLASMA_request*,     request);
            }
        }
    }

    return PLASMA_SUCCESS;
}
