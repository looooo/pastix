/**
 *
 * @file async.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#include "common.h"

#include <stdlib.h>

/***************************************************************************//**
 *  Register an exception.
 **/
int plasma_request_fail(PLASMA_sequence *sequence, PLASMA_request *request, int status)
{
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return status;
}

/***************************************************************************//**
 *  Create a sequence
 **/
int plasma_sequence_create(plasma_context_t *plasma, PLASMA_sequence **sequence)
{
    if ((*sequence = malloc(sizeof(PLASMA_sequence))) == NULL) {
        plasma_error("PLASMA_Sequence_Create", "malloc() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    if(((*sequence)->quark_sequence = QUARK_Sequence_Create(plasma->quark)) == NULL){
        plasma_error("PLASMA_Sequence_Create", "QUARK_Sequence_Create() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    (*sequence)->status = PLASMA_SUCCESS;
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *  Destroy a sequence
 **/
int plasma_sequence_destroy(plasma_context_t *plasma, PLASMA_sequence *sequence)
{
    QUARK_Sequence_Destroy(plasma->quark, sequence->quark_sequence);
    free(sequence);
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *  Wait for the completion of a sequence
 **/
int plasma_sequence_wait(plasma_context_t *plasma, PLASMA_sequence *sequence)
{
    QUARK_Sequence_Wait(plasma->quark, sequence->quark_sequence);
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *  Terminate a sequence
 **/
void plasma_sequence_flush(Quark *quark, PLASMA_sequence *sequence, PLASMA_request *request, int status)
{
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    QUARK_Sequence_Cancel(quark, sequence->quark_sequence);
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Sequence_Create - Create a squence.
 *
 *******************************************************************************
 *
 * @param[out] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Sequence_Create(PLASMA_sequence **sequence)
{
    plasma_context_t *plasma;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Sequence_Create", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    status = plasma_sequence_create(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Sequence_Destroy - Destroy a sequence.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Sequence_Destroy(PLASMA_sequence *sequence)
{
    plasma_context_t *plasma;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Sequence_Destroy", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_Sequence_Destroy", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    status = plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Sequence_Wait - Wait for the completion of a sequence.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Sequence_Wait(PLASMA_sequence *sequence)
{
    plasma_context_t *plasma;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Sequence_Wait", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_Sequence_Wait", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    status = plasma_sequence_wait(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Sequence_Flush - Terminate a sequence.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 * @param[in] request
 *          The flush request.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Sequence_Flush(PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Sequence_Flush", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_Sequence_Flush", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    plasma_sequence_flush(plasma->quark, sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);
    return PLASMA_SUCCESS;
}
