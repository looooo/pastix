/**
 *
 * @file context.c
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
#include <stdlib.h>
#if defined( _WIN32 ) || defined( _WIN64 )
#include "plasmawinthread.h"
#else
#include <pthread.h>
#endif

#include "common.h"
#include "auxiliary.h"
#include "context.h"

/***************************************************************************//**
 *  Global data
 **/
/* master threads context lookup table */
plasma_context_map_t context_map[CONTEXTS_MAX];
/* context lookup table access lock */
pthread_mutex_t context_map_lock = PTHREAD_MUTEX_INITIALIZER;

/***************************************************************************//**
 *  Create new context
 **/
plasma_context_t *plasma_context_create()
{
    plasma_context_t *plasma;

    plasma = (plasma_context_t*)malloc(sizeof(plasma_context_t));
    if (plasma == NULL) {
        plasma_fatal_error("plasma_context_create", "malloc() failed");
        return NULL;
    }

    pthread_mutex_init(&plasma->action_mutex, NULL);
    pthread_mutex_init(&context_map_lock, NULL);
    pthread_cond_init(&plasma->action_condt, NULL);
    plasma->action = PLASMA_ACT_STAND_BY;
    plasma->parallel_func_ptr = NULL;

    /* These initializations are just in case the user
       disables autotuning and does not set nb and ib */
    plasma->nb = 128;
    plasma->ib = 32;
    plasma->nbnbsize = 16384;
    plasma->ibnbsize = 4096;
    plasma->rhblock = 4;

    plasma->errors_enabled     = PLASMA_TRUE;
    plasma->warnings_enabled   = PLASMA_TRUE;
    plasma->autotuning_enabled = PLASMA_TRUE;
    plasma->dynamic_section    = PLASMA_FALSE;

    plasma->scheduling  = PLASMA_STATIC_SCHEDULING;
    plasma->householder = PLASMA_FLAT_HOUSEHOLDER;
    plasma->translation = PLASMA_OUTOFPLACE;

    return plasma;
}

/***************************************************************************//**
 *  Insert a (context, thread_id) tuple in the context map
 **/
int plasma_context_insert(plasma_context_t *context, pthread_t thread_id)
{
    int i;

    // acquire the access lock
    pthread_mutex_lock(&context_map_lock);
    // For each entry
    for (i = 0; i < CONTEXTS_MAX; i++) {
        // If not occupied
        if (context_map[i].context == NULL) {
            // Insert new context, release lock, return success
            context_map[i].context = context;
            context_map[i].thread_id = thread_id;
            pthread_mutex_unlock(&context_map_lock);
            return PLASMA_SUCCESS;
        }
    }
    // No empty entry found - release lock, return error
    pthread_mutex_unlock(&context_map_lock);
    plasma_fatal_error("plasma_context_insert", "too many threads");
    return PLASMA_ERR_INTERNAL_LIMIT;
}

/***************************************************************************//**
 *  Remove a (context, thread_id) tuple from the context map
 **/
int plasma_context_remove(plasma_context_t *context, pthread_t thread_id)
{
    int i;

    // acquire the access lock
    pthread_mutex_lock(&context_map_lock);
    // For each entry
    for (i = 0; i < CONTEXTS_MAX; i++) {
        // If id matches
        if (pthread_equal(context_map[i].thread_id, thread_id)) {
            if (context_map[i].context == context) {
                // Free the context, mark entry as empty, release lock, return success
                free(context_map[i].context);
                context_map[i].context = NULL;
                pthread_mutex_unlock(&context_map_lock);
                return PLASMA_SUCCESS;
            }
            else {
                pthread_mutex_unlock(&context_map_lock);
                plasma_fatal_error("plasma_context_remove", "context does not match thread");
                return PLASMA_ERR_UNEXPECTED;
            }
        }
    }
    // No matching id found - release lock, return error
    pthread_mutex_unlock(&context_map_lock);
    plasma_fatal_error("plasma_context_remove", "thread not found");
    return PLASMA_ERR_NOT_FOUND;
}

/***************************************************************************//**
 *  Return context for a thread
 **/
plasma_context_t *plasma_context_self()
{
    int i;

    // For each entry
    for (i = 0; i < CONTEXTS_MAX; i++) {
        // If id matches
        if (pthread_equal(context_map[i].thread_id, pthread_self())) {
            return context_map[i].context;
        }
    }
    return NULL;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Enable - Enable PLASMA feature.
 *
 *******************************************************************************
 *
 * @param[in] lever
 *          Feature to be enabled:
 *          @arg PLASMA_WARNINGS   printing of warning messages,
 *          @arg PLASMA_ERRORS     printing of error messages,
 *          @arg PLASMA_AUTOTUNING autotuning for tile size and inner block size,
 *          @arg PLASMA_DAG        Quark automatic DAG generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Enable(PLASMA_enum lever)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Enable", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    switch (lever)
    {
        case PLASMA_WARNINGS:
            plasma->warnings_enabled = PLASMA_TRUE;
            break;
        case PLASMA_ERRORS:
            plasma->errors_enabled = PLASMA_TRUE;
            break;
        case PLASMA_AUTOTUNING:
            plasma->autotuning_enabled = PLASMA_TRUE;
            break;
        case PLASMA_DAG:
            QUARK_Barrier(plasma->quark);
            QUARK_DOT_DAG_Enable( plasma->quark, 1);
            break;
        default:
            plasma_error("PLASMA_Enable", "illegal parameter value");
            return PLASMA_ERR_ILLEGAL_VALUE;
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Disable - Disable PLASMA feature.
 *
 *******************************************************************************
 *
 * @param[in] lever
 *          Feature to be disabled:
 *          @arg PLASMA_WARNINGS   printing of warning messages,
 *          @arg PLASMA_ERRORS     printing of error messages,
 *          @arg PLASMA_AUTOTUNING autotuning for tile size and inner block size,
 *          @arg PLASMA_DAG        Quark automatic DAG generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Disable(PLASMA_enum lever)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Disable", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    switch (lever)
    {
        case PLASMA_WARNINGS:
            plasma->warnings_enabled = PLASMA_FALSE;
            break;
        case PLASMA_ERRORS:
            plasma->errors_enabled = PLASMA_FALSE;
            break;
        case PLASMA_AUTOTUNING:
            plasma->autotuning_enabled = PLASMA_FALSE;
            break;
        case PLASMA_DAG:
            QUARK_Barrier(plasma->quark);
            QUARK_DOT_DAG_Enable( plasma->quark, 0);
            break;
        default:
            plasma_error("PLASMA_Disable", "illegal parameter value");
            return PLASMA_ERR_ILLEGAL_VALUE;
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Set - Set PLASMA parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg PLASMA_TILE_SIZE:        size matrix tile,
 *          @arg PLASMA_INNER_BLOCK_SIZE: size of tile inner block,
 *          @arg PLASMA_SCHEDULING_MODE:  type of scheduling
 *                  (PLASMA_STATIC_SCHEDULING or PLASMA_DYNAMIC_SCHEDULING).
 *          @arg PLASMA_INNER_BANDWIDTH:  matrix bandwidth size
 *
 * @param[in] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Set(PLASMA_enum param, int value)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA_Set", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case PLASMA_TILE_SIZE:
            if (value <= 0) {
                plasma_error("PLASMA_Set", "negative tile size");
                return PLASMA_ERR_ILLEGAL_VALUE;
            }
            plasma->nb = value;
            /* Calculate A, B tile size and round up to cache line size */
            /* round up for the smallest type (float) - will hold for all */
            plasma->nbnbsize = plasma->nb * plasma->nb; // * sizeof(float);
            plasma->ibnbsize = plasma->ib * plasma->nb; // * sizeof(float);
//          plasma->nbnbsize = roundup(plasma->nbnbsize, CACHE_LINE_SIZE);
//          plasma->nbnbsize /= sizeof(float);
            if ( plasma->autotuning_enabled ) {
                plasma->autotuning_enabled = PLASMA_FALSE;
                plasma_warning("PLASMA_Set", "autotuning has been automatically disabled\n");
            }
            /* Limit ib to nb */
            plasma->ib = min( plasma->nb, plasma->ib );
            break;
        case PLASMA_INNER_BLOCK_SIZE:
            if (value <= 0) {
                plasma_error("PLASMA_Set", "negative inner block size");
                return PLASMA_ERR_ILLEGAL_VALUE;
            }
            if (value > plasma->nb) {
                plasma_error("PLASMA_Set", "inner block larger than tile");
                return PLASMA_ERR_ILLEGAL_VALUE;
            }
            /* if (plasma->nb % value != 0) { */
            /*     plasma_error("PLASMA_Set", "inner block does not divide tile"); */
            /*     return PLASMA_ERR_ILLEGAL_VALUE; */
            /* } */
            plasma->ib = value;
            /* Calculate T, L tile size and round up to cache line size */
            /* round up for the smallest type (float) - will hold for all */
            plasma->ibnbsize = plasma->ib * plasma->nb; // * sizeof(float);
//          plasma->ibnbsize = roundup(plasma->ibnbsize, CACHE_LINE_SIZE);
//          plasma->ibnbsize /= sizeof(float);

            if ( plasma->autotuning_enabled ) {
                plasma->autotuning_enabled = PLASMA_FALSE;
                plasma_warning("PLASMA_Set", "autotuning has been automatically disabled\n");
            }
            break;
        case PLASMA_SCHEDULING_MODE:
            if (value != PLASMA_STATIC_SCHEDULING && value != PLASMA_DYNAMIC_SCHEDULING) {
                plasma_error("PLASMA_Set", "illegal value of PLASMA_SCHEDULING_MODE");
                return PLASMA_ERR_ILLEGAL_VALUE;
            }
            plasma->scheduling = value;
            break;
        case PLASMA_HOUSEHOLDER_MODE:
            if (value != PLASMA_FLAT_HOUSEHOLDER && value != PLASMA_TREE_HOUSEHOLDER) {
                plasma_error("PLASMA_Set", "illegal value of PLASMA_HOUSEHOLDER_MODE");
                return PLASMA_ERR_ILLEGAL_VALUE;
            }
            plasma->householder = value;
            break;
        case PLASMA_HOUSEHOLDER_SIZE:
            if (value <= 0) {
                plasma_error("PLASMA_Set", "negative householder size");
                return PLASMA_ERR_ILLEGAL_VALUE;
            }
            plasma->rhblock = value;
            break;
        case PLASMA_TRANSLATION_MODE:
            if (value != PLASMA_INPLACE && value != PLASMA_OUTOFPLACE) {
                plasma_error("PLASMA_Set", "illegal value of PLASMA_TRANSLATION_MODE");
                return PLASMA_ERR_ILLEGAL_VALUE;
            }
            plasma->translation = value;
            break;
        default:
            plasma_error("PLASMA_Set", "unknown parameter");
            return PLASMA_ERR_ILLEGAL_VALUE;
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Get - Get value of PLASMA parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg PLASMA_TILE_SIZE:        size matrix tile,
 *          @arg PLASMA_INNER_BLOCK_SIZE: size of tile inner block,
 *          @arg PLASMA_SCHEDULING:       type of scheduling
 *                  (PLASMA_STATIC_SCHEDULING or PLASMA_DYNAMIC_SCHEDULING).
 *
 * @param[out] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Get(PLASMA_enum param, int *value)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Get", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case PLASMA_TILE_SIZE:
            *value = plasma->nb;
            return PLASMA_SUCCESS;
        case PLASMA_INNER_BLOCK_SIZE:
            *value = plasma->ib;
            return PLASMA_SUCCESS;
        case PLASMA_SCHEDULING_MODE:
            *value = plasma->scheduling;
            return PLASMA_SUCCESS;
        case PLASMA_HOUSEHOLDER_MODE:
            *value = plasma->householder;
            return PLASMA_SUCCESS;
        case PLASMA_HOUSEHOLDER_SIZE:
            *value = plasma->rhblock;
            return PLASMA_SUCCESS;
        case PLASMA_TRANSLATION_MODE:
            *value = plasma->translation;
            return PLASMA_SUCCESS;
        default:
            plasma_error("PLASMA_Get", "unknown parameter");
            return PLASMA_ERR_ILLEGAL_VALUE;
    }
}
